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
	int	    nPoints;



double integral(double x, double y, double z,std::vector<Point> deposits){
	//absorbtionTerm will compute the absorbed ammount of power from the laser
	//ContributionSum will sum up contributions from all deposits
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
    return contributionSum*absorbtionTerm*dv/(4*pi*waterConductivity); 
}


int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    = stold(argv[2])  * pow(10,-6); 
	stepSize  = bounds/(double)(300);		  //Step size, based off of bounds parameter
	nDeposits = stof(argv[1]);				  //number of deposits to initialize
	//size of simulation box
	lambda	       = stold(argv[3])  * pow(10,-9); //Spatial periodicity
    dv       	   = stepSize*stepSize*stepSize;  //volume element for integral
	
	int nParticles =1;

	Point centerOfParticle1 = {0.0,0.0,0.0}; 
	Point centerOfParticle2 = {-1.05*particleRadius,0.0,0.0}; 
	Particle particle1(centerOfParticle1,particleRadius);
	Particle particle2(centerOfParticle2,particleRadius);
	vector<Particle> particles = {particle1,particle2};
     ///////////GENERATE DEPOSITS//////////////////////////////////////////////////////////////////
	 for(int i = 0; i<nParticles; i++ ){
		particles[i].generateDeposits(nDeposits);
	 }

	//////////////////////////////////////////////////////////////////////////////////////////////
	///////////INITIALIZE LINSPACE VECTORS////////////////////////////////////////////////////////
	vector<double> linspace = arange(-bounds,bounds,stepSize);
	vector<double> y = {0.0};
    vector<double> z = linspace;
	vector<double> x = linspace;

   
 	cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
 	//////////////////////////////////////////////////////////////////////////////////////////////
 	////////////////////////  INTEGRAL ///////////////////////////////////////////////////////////

	nPoints = x.size();
	vector<vector<vector<double>>> xGrad(nPoints, vector<vector<double>>(nPoints, vector<double>(nPoints)));
    vector<vector<vector<double>>> zGrad(nPoints, vector<vector<double>>(nPoints, vector<double>(nPoints)));
      
    const int totalIterations = nPoints*nPoints*y.size();
	size_t currentIteration = 0;

	/////////////////////////////////
	//Compute gradient in X-direction
	/////////////////////////////////
	double dl = stepSize*stepSize;
for(int n = 0; n < nParticles;n++){

	double thickness = pow(100*stepSize,2);
	#pragma omp parallel for
    for (size_t i = 1; i < nPoints-1; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 1; k<nPoints-1; k++){		
    			//Check if outside particle:
				Point point = {x[i],y[j],z[k]};
				double d = particles[n].getRadialDistance(point);
				if (d > pow(particles[n].radius,2)){
					
					//double back   			  = integral(x[i-1],y[j],z[k],particle.deposits);
					//double forward			  = integral(x[i+1],y[j],z[k],particle.deposits);

					double back   			  = integral(x[i]-dl,y[j],z[k],particles[n].deposits);
					double forward			  = integral(x[i]+dl,y[j],z[k],particles[n].deposits);

					double central_difference = (forward - back)/(2*stepSize);

					double perpendicular      = central_difference*x[i]*x[i]/d;

					double tangential         = central_difference - perpendicular;
					//xGrad[i][j][k] 	      += perpendicular*25/1000;	
					xGrad[i][j][k] 			  -= tangential*25/1000;		//Converts to Kelvin/micrometer 
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
    for (size_t i = 1; i < nPoints-1; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 1; k<nPoints-1; k++){		
    			//Check if outside particle:
				Point point = {x[i],y[j],z[k]};
				double d = particles[n].getRadialDistance(point);
				if (d > pow(particles[n].radius,2)){
					
					double back   			  = integral(x[i],y[j],z[k]-dl,particles[n].deposits);
					double forward 			  = integral(x[i],y[j],z[k]+dl,particles[n].deposits);
					double central_difference = (forward - back)/(2*stepSize);

					double perpendicular      = central_difference*z[k]*z[k]/d;
					double tangential         = central_difference - perpendicular;

				//	zGrad[i][j][k] 		  	  += perpendicular*25/1000;	
					zGrad[i][j][k] 			  += tangential*25/1000;	//Converts to Kelvin/micrometer 	
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
	particles[n].writeDepositToCSV();
}

    //////////////////////////////////////////////////////////////////
    //////////////////////WRITE TO FILE///////////////////////////////
	cout<<"Simulation finished, writing to csv..."<<endl;
	writeGradToCSV(x,y,z,xGrad,zGrad);
	
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