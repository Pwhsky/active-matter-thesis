#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
#include "functions.h"


using namespace std;
 	double  bounds;
	double  lambda;
	double  stepSize;
	double  dv;
	int 	nDeposits;	
	int	nSteps;



double integral(double x, double y, double z,std::vector<Point> deposits){
	//absorbtionTerm will compute the absorbed ammount of power from the laser
	//ContributionSum will sum up contributions
	//Finally, the contributionSum is scaled with volume element dv and divided with constants												
	double laserPower	           = I0 + I0*cos(twoPi*(x)/lambda);	

	double absorbtionTerm          = laserPower*depositArea/(volumePerDeposit);
	double contributionSum 		   = 0.0;
		//Since the values scale with the inverse square distance, this is the only heavy computation.
    	for (size_t i = 0; i < deposits.size(); i++){
    		double inverse_squareroot_distance = 1.0/sqrt(pow(x-deposits[i].x,2)+
												 pow(y-deposits[i].y,2)+
												 pow(z-deposits[i].z,2));
			contributionSum +=  inverse_squareroot_distance;
		}
    return contributionSum*dv*absorbtionTerm/(4*pi*waterConductivity); 
}


int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    = stold(argv[2])  * pow(10,-6); 
	stepSize  = bounds/(double)(300);		  //Step size, based off of bounds parameter
	nDeposits = stof(argv[1]);				  //number of deposits to initialize
	//size of simulation box
	lambda	  = stold(argv[3])  * pow(10,-9); //Spatial periodicity
    dv        = stepSize*stepSize*stepSize;  //volume element for integral
    
	Point centerOfParticle = {0.0,0.0,particleRadius}; 
	Particle particle(centerOfParticle,particleRadius);
     ///////////GENERATE DEPOSITS//////////////////////////////////////////////////////////////////
	particle.generateDeposits(particle.deposits,nDeposits);

	//////////////////////////////////////////////////////////////////////////////////////////////
	///////////INITIALIZE LINSPACE VECTORS////////////////////////////////////////////////////////
	vector<double> linspace = arange(-bounds,bounds,stepSize);
	vector<double> y = {0.0};
    vector<double> z = linspace;
	vector<double> x = linspace;

   
 	cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
 	//////////////////////////////////////////////////////////////////////////////////////////////
 	////////////////////////  INTEGRAL ///////////////////////////////////////////////////////////

	nSteps = x.size();
	vector<vector<vector<double>>> xGrad(nSteps, vector<vector<double>>(nSteps, vector<double>(nSteps)));
    vector<vector<vector<double>>> zGrad(nSteps, vector<vector<double>>(nSteps, vector<double>(nSteps)));
      
    const int totalIterations = nSteps*nSteps*y.size();
	size_t currentIteration = 0;

	/////////////////////////////////
	//Compute gradient in X-direction
	/////////////////////////////////
	#pragma omp parallel for
    for (size_t i = 1; i < nSteps-1; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 1; k<nSteps-1; k++){		
    			//Check if outside particle:
				Point point = {x[i],y[j],z[k]};
				if (particle.isOutside(point)){

					double back   			  = integral(x[i-1],y[j],z[k],particle.deposits);
					double forward			  = integral(x[i+1],y[j],z[k],particle.deposits);
					double central_difference = (forward - back)/(2*stepSize);
					xGrad[i][j][k] 			  = central_difference*25/1000;		//Converts to Kelvin/micrometer 
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
	/////////////////////////////////
	//Compute gradient in Z-direction
	/////////////////////////////////
	#pragma omp parallel for
    for (size_t i = 1; i < nSteps-1; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 1; k<nSteps-1; k++){		
    			//Check if outside particle:
				Point point = {x[i],y[j],z[k]};

				if (particle.isOutside(point)){
					double back   			  = integral(x[i],y[j],z[k-1],particle.deposits);
					double forward 			  = integral(x[i],y[j],z[k+1],particle.deposits);
					double central_difference = (forward - back)/(2*stepSize);
					zGrad[i][j][k] 			  = central_difference*25/1000;	//Converts to Kelvin/micrometer 	
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

    //////////////////////////////////////////////////////////////////
    //////////////////////WRITE TO FILE///////////////////////////////
	cout<<"Simulation finished, writing to csv..."<<endl;
	writeGradToCSV(x,y,z,xGrad,zGrad);
	writeDepositToCSV(particle.deposits);
	
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