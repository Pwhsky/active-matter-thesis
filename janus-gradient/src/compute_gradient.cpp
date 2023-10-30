#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
#include "functions.h"
#include <cstring>

using namespace std;
 	double  bounds;
	double  lambda;
	double  stepSize;
	double  dv;
	int 	nDeposits;	
	int	nSteps;

inline  double integral(double x, double y, double z,vector<Point> deposits){
	//Pre-compute the terms, laserTerm will give power for a given displacement from the center
	//absorbtionTerm will compute the absorbed ammount of power from laserTerm.
	//ContributionSum will sum up contributions
	//Finally, the contributionSum is scaled with volume element dv and divided with constants
	//		
	double laserTerm 	        = I0 + I0*cos(twoPi*(x)/lambda);
	double absorbtionTerm       = laserTerm*depositArea/(volumePerDeposit);
	double contributionSum 		= 0.0;
	double q                	= 0.0;
    for (size_t i = 0; i < deposits.size(); i++){
    	double inv_sqrt_distance1 = 1.0/sqrt((x-deposits[i].x)*(x-deposits[i].x) +  (y-deposits[i].y)*(y-deposits[i].y) +(z-deposits[i].z)*(z-deposits[i].z));
  		contributionSum +=  inv_sqrt_distance1;
	}
    return contributionSum*dv*absorbtionTerm/(4*pi*waterConductivity);
}


int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    = stold(argv[3])  * pow(10,-6); 
	stepSize  = bounds/(stof(argv[1]));		  //Step size, based off of resolution parameter
	nDeposits = stof(argv[2]);				  //number of deposits to initialize
	//size of simulation box
	lambda	  = stold(argv[4])  * pow(10,-9); //Spatial periodicity
    dv        = stepSize*stepSize*stepSize;  //volume element for integral
        
     ///////////GENERATE DEPOSITS//////////////////////////////////////////////////////////////////
	vector<Point> deposits;
	vector<Point> boundaryPoints;
	generateDeposits(deposits,nDeposits);
	//////////////////////////////////////////////////////////////////////////////////////////////
	///////////INITIALIZE LINSPACE VECTORS////////////////////////////////////////////////////////
	vector< double> linspace;
	vector<double> theta;
	double coordinate;
	vector<double> y;
	for (coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
          	 linspace.push_back(coordinate);
    }

	for (coordinate = 0; coordinate <= twoPi; coordinate += twoPi/500) {
			Point point;
			point.x = cos(coordinate)*particleRadius;
			point.z = sin(coordinate)*particleRadius;
			point.y = 0.0;
			boundaryPoints.push_back(point);
    }

    vector<double> z = linspace;
	vector<double> x = linspace;
	//y= linspace;
	y = {0.0}; 
   
 	cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
 	//////////////////////////////////////////////////////////////////////////////////////////////
 	////////////////////////  INTEGRAL ///////////////////////////////////////////////////////////

	nSteps = x.size();
	vector<vector<vector<double>>> xGrad(nSteps, vector<vector<double>>(nSteps, vector<double>(nSteps)));
    vector<vector<vector<double>>> zGrad(nSteps, vector<vector<double>>(nSteps, vector<double>(nSteps)));
      
    const int totalIterations = nSteps*nSteps*y.size();
	size_t currentIteration = 0;

	/*
	//Goal is to compute Radial and Tangental flow 
 	#pragma omp parallel for
    for (size_t i = 0; i < nSteps; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 0; k<nSteps; k++){		
    			//Check if outside particle:

				double contribution = 0.0;
				double greensNorm = pow(x[i],2)+pow(y[j],2) + pow(z[k],2);
				if (greensNorm > particleRadiusSquared){
				for (size_t ix = 0; ix < boundaryPoints.size(); ix++ ){
					Point bdp = boundaryPoints[ix];

					double norm = bdp.x*bdp.x + bdp.y*bdp.y + bdp.z*bdp.z;
					double distance = pow(x[i]-bdp.x,2) + pow(y[j]-bdp.y,2) + pow(z[k]-bdp.z,2);
			
					double back    = integral(bdp.x,bdp.y,(bdp.z-stepSize/2),deposits);
					double forward = integral(bdp.x,bdp.y,(bdp.z+stepSize/2),deposits);
					double difference = (forward - back)/(stepSize);
					double parallel   = difference - difference*(pow(bdp.z,2) + bdp.z)/norm; 
					contribution -= parallel/sqrt(distance);
				}
					zGrad[i][j][k] = contribution*stepSize*stepSize/8;
				}
			
				currentIteration++;
         		// Calculate progress percentage so that the user has something to look at
          		if(currentIteration % 500 == 0) {
            	  	float progress = round(static_cast<float>(currentIteration) / totalIterations * 100.0);
            		 // Print progress bar
            		 #pragma omp critical
            		{
            		 		cout << "Progress: "<< progress << "% ("<< currentIteration<< "/" << totalIterations << ")\r";
            	  	 		cout.flush();
            	 	}
            	}
			}
    	}
    }
	 	#pragma omp parallel for
    for (size_t i = 0; i < nSteps; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 0; k<nSteps; k++){		
    			//Check if outside particle:

				double contribution = 0.0;
				double greensNorm = pow(x[i],2)+pow(y[j],2) + pow(z[k],2);
				if (greensNorm > particleRadiusSquared){
				for (size_t ix = 0; ix < boundaryPoints.size(); ix++ ){
					Point bdp = boundaryPoints[ix];

					double norm = bdp.x*bdp.x + bdp.y*bdp.y + bdp.z*bdp.z;
					double distance = pow(x[i]-bdp.x,2) + pow(y[j]-bdp.y,2) + pow(z[k]-bdp.z,2);
			
					double back    = integral(bdp.x-stepSize/2,bdp.y,(bdp.z),deposits);
					double forward = integral(bdp.x+stepSize/2,bdp.y,(bdp.z),deposits);
					double difference = (forward - back)/(stepSize);
					double parallel   = difference - difference*(pow(bdp.x,2) + bdp.x)/norm; 
					contribution -= parallel/sqrt(distance);
				}
					xGrad[i][j][k] = contribution*stepSize*stepSize/8;
				}
			
				currentIteration++;
         		// Calculate progress percentage so that the user has something to look at
          		if(currentIteration % 500 == 0) {
            	  	float progress = round(static_cast<float>(currentIteration) / totalIterations * 100.0);
            		 // Print progress bar
            		 #pragma omp critical
            		{
            		 		cout << "Progress: "<< progress << "% ("<< currentIteration<< "/" << totalIterations << ")\r";
            	  	 		cout.flush();
            	 	}
            	}
			}
    	}
    } */
//Goal is to compute Radial and Tangental flow 
	#pragma omp parallel for
    for (size_t i = 1; i < nSteps-1; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 1; k<nSteps-1; k++){		
    			//Check if outside particle:
				if (x[i]*x[i] + y[j]*y[j] + z[k]*z[k] > particleRadiusSquared)	{
					double back    = integral(x[i-1],y[j],z[k],deposits);
					double forward = integral(x[i+1],y[j],z[k],deposits);
					double difference = (forward - back)/(2*stepSize);
					xGrad[i][j][k] = difference*25/1000;		
      			}
				currentIteration++;
         		// Calculate progress percentage so that the user has something to look at
          		if(currentIteration % 500 == 0) {
            	  	float progress = round(static_cast<float>(currentIteration) / totalIterations * 100.0);
            		 // Print progress bar
            		 #pragma omp critical
            		{
            		 		cout << "Progress: "<< progress << "% ("<< currentIteration<< "/" << totalIterations << ")\r";
            	  	 		cout.flush();
            	 	}
            	}
			}
    	}
    }
	#pragma omp parallel for
    for (size_t i = 1; i < nSteps-1; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 1; k<nSteps-1; k++){		
    			//Check if outside particle:
				if (x[i]*x[i] + y[j]*y[j] + z[k]*z[k] > particleRadiusSquared)	{
					double back    = integral(x[i],y[j],z[k-1],deposits);
					double forward = integral(x[i],y[j],z[k+1],deposits);
					double difference = (forward - back)/(2*stepSize);
					zGrad[i][j][k] = difference*25/1000;		
      			}
				currentIteration++;
         		// Calculate progress percentage so that the user has something to look at
          		if(currentIteration % 500 == 0) {
            	  	float progress = round(static_cast<float>(currentIteration) / totalIterations * 100.0);
            		 // Print progress bar
            		 #pragma omp critical
            		{
            		 		cout << "Progress: "<< progress << "% ("<< currentIteration<< "/" << totalIterations << ")\r";
            	  	 		cout.flush();
            	 	}
            	}
			}
    	}
    }


	//Find out a way to convert x and y to theta and r
	cout<<"Finished computing boundary, starting flow field computation... \n";

    	//computeGradientZ(field, x,  y, z);

	/*
 	#pragma omp parallel for
    for (size_t i = 1; i < nSteps-1; i++){
    	for(size_t j = 0; j < y.size(); j++){
    		for(size_t k = 1; k<nSteps-1; k++){
    		
    			//Check if outside particle:
				if (x[i]*x[i] + y[j]*y[j] + z[k]*z[k] > particleRadiusSquared)	{
					double back       = integral(x[i],y[j],z[k+1],deposits);
					double forward    = integral(x[i],y[j],z[k-1],deposits);
					double difference = (forward - back)/(2*stepSize);
					field[i][j][k]    += pow((difference*25/1000),2); //normalize to kelvin/micrometer
					field[i][j][k]    = sqrt(field[i][j][k]);
       			}
				currentIteration++;
			    	// Calculate progress percentage so that the user has something to look at
              	if(currentIteration % 500 == 0) {
             		float progress = round(static_cast<float>(currentIteration) / totalIterations * 100.0);
               		// Print progress bar
               		#pragma omp critical
               		{
                		cout << "Progress: "<< progress << "% ("<< currentIteration<< "/" << totalIterations << ")\r";
                	  	cout.flush();
                	}
              	}
			}
    	}
    }
    		*/
    	//////////////////////////////////////////////////////////////////
    	//////////////////////WRITE TO FILE///////////////////////////////
	cout<<"Simulation finished, writing to csv..."<<endl;
	writeGradToCSV(x,y,z,xGrad,zGrad);
	writeDepositToCSV(deposits);	
	
	
	
	//////////////////////////////////////////////////////////////////
	///////////////////COMPUTE ELAPSED TIME///////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Program completed after: " << elapsed_seconds << " seconds" << std::endl;
	//////////////////////////////////////////////////////////////////	
	//////////////////END PROGRAM/////////////////////////////////////
	return 0;
}	



