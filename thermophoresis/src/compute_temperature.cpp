#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
#include "particle.h"



using namespace std;
 	double bounds, lambda, stepSize, dv;
	int    nDeposits;	

int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();

	//Parse input arguments:
	bounds    = stold(argv[2])  * pow(10,-6); //size of simulation box
	stepSize  = bounds/(double)(300);       //step size based off of bounds parameter
	nDeposits = stof(argv[1]);				  //number of deposits to initialize
	lambda	  = stold(argv[3])  * pow(10,-9); //Spatial periodicity of laser 
    dv        = stepSize*stepSize*stepSize;   //volume element for integral
	

	vector<Particle> particles = initializeParticles();
	int nParticles = particles.size();
	for(int i = 0; i<nParticles; i++ ) particles[i].generateDeposits(nDeposits);

	///////////INITIALIZE LINSPACE VECTORS////////////////////////////////////////////////////////
	 vector<double> linspace = arange(-1*bounds,bounds,stepSize);

     const vector<double> z = linspace;
	 const vector<double> x = linspace;
 	 vector<double>       y = {0.0};
	 //vector<double>	  y = linspace; 
 	 cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
	 
 	//////////////////////////////////////////////////////////////////////////////////////////////
 	////////////////////////  INTEGRAL ///////////////////////////////////////////////////////////
	int nSteps = x.size();
	//Initialize empty 3D space of points which will contain the integral values.
	 vector<vector<vector<double>>> temperature(nSteps, vector<vector<double>>(nSteps, vector<double>(nSteps)));  
	const int totalIterations = nSteps*nSteps*y.size();
   	size_t currentIteration = 0;

	//The 3 nested for-loops will call integral() for all points in 3D space.
	for(int n = 0;n<nParticles;n++){
 	#pragma omp parallel for
    	for (size_t i = 0; i < nSteps; i++){
    		for(size_t j = 0; j<y.size(); j++){
    			for(size_t k = 0; k<nSteps; k++){
					Point point = {x[i],y[j],z[k]};
					double d = particles[n].getRadialDistance(point);
    				//Check if outside particle:
					if (d > pow(particles[n].radius,2)){
						temperature[i][j][k] += integral(x[i],y[j],z[k],particles[n].deposits,lambda,dv);	
       				}

					// Display progress bar so that the user has something to look at
					currentIteration++;
              		if(currentIteration % 500 == 0) {
             			float progress = round(static_cast<float>(currentIteration) / totalIterations * 100.0);
               			#pragma omp critical
               			{
								//Print progress bar
                				cout << "Progress: "<< progress << "% ("<< currentIteration<< "/" << totalIterations << ")\r";
                	  	 		cout.flush();
                		}
              		}
				}
    		}
    	}
		particles[n].writeDepositToCSV();
	}

    //////////////////////WRITE TO FILE///////////////////////////////
		cout<<"Simulation finished, writing to csv..."<<endl;
		writeFieldToCSV(x,y,z,temperature);	
	///////////////////COMPUTE ELAPSED TIME///////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Program completed after: " << elapsed_seconds << " seconds" << std::endl;
	//////////////////END PROGRAM/////////////////////////////////////
	return 0;
}	


