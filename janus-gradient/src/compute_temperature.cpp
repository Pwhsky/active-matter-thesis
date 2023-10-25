#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
#include "functions.h"
#include <cstring>


using namespace std;
 	double bounds;
	double lambda;
	double stepSize;
	double dv;
	int 	    nDeposits;	

inline  double integral(double x, double y, double z,vector<Point> deposits){

	//Pre-compute the terms, laserTerm will give power for a given displacement from the center
	//absorbtionTerm will compute the absorbed ammount of power from laserTerm.
	//ContributionSum will sum up contributions
	//Finally, the contributionSum is scaled with volume element dv and divided with constants
	//																	to get physical results.

	double laserTerm 	           = I0 + I0*cos(twoPi*(x)/lambda);	
	double absorbtionTerm          = laserTerm*depositArea/(volumePerDeposit);
	double contributionSum 		   = 0.0;
		//Since the values scale with the inverse square distance, this is the only heavy computation.
    	for (size_t i = 0; i < deposits.size(); i++){
    		double inv_sqrt_distance1 = 1.0/sqrt((x-deposits[i].x)*(x-deposits[i].x) + (y-deposits[i].y)*(y-deposits[i].y) + z-deposits[i].z)*(z-deposits[i].z));
			contributionSum +=  inv_sqrt_distance1;
		}
    return contributionSum*dv*absorbtionTerm/(4*pi*waterConductivity);
}




int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();

	//Parse input arguments:
	stepSize  = bounds/(stof(argv[1]));     //step size based off of resolution parameter, 200 = 25 nm step size
	nDeposits = stof(argv[2]);				//number of deposits to initialize
	bounds    = stold(argv[3])  * pow(10,-6); //size of simulation box		 in micrometers
	lambda	  = stold(argv[4])  * pow(10,-9); //Spatial periodicity of laser in nanometers
    dv        = stepSize*stepSize*stepSize; //volume element for integral
	
	
	///////////GENERATE DEPOSITS//////////////////////////////////////////////////////////////////
	vector<Point> deposits;
	generateDeposits(deposits,nDeposits);
	//generateConfiguration(deposits,nDeposits);
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	///////////INITIALIZE LINSPACE VECTORS////////////////////////////////////////////////////////
	 vector< double> spaceVector;
	 for (double coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
          	 spaceVector.push_back(coordinate);
      	 }
     	 const vector<double> z = spaceVector;
	 const vector<double> x = spaceVector;
	 //vector<double>	y = spaceVector;
	 vector<double>       y = {0.0};
    
 	 cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
 	//////////////////////////////////////////////////////////////////////////////////////////////
 	////////////////////////  INTEGRAL ///////////////////////////////////////////////////////////
	int nSteps = x.size();
	//Initialize empty 3D space of points which will contain the integral values.
	 vector<vector<vector<double>>> field(nSteps, vector<vector<double>>(nSteps, vector<double>(nSteps)));  
	const int totalIterations = nSteps*nSteps;
   	size_t currentIteration = 0;
 	#pragma omp parallel for
		//The 3 for-loops will check all points in 3D space and compute the integral
    	for (size_t i = 0; i < nSteps; i++){
    		for(size_t j = 0; j<y.size(); j++){
    			for(size_t k = 0; k<nSteps; k++){
    			
    				//Check if outside particle:
					if (x[i]*x[i] + y[j]*y[j] + z[k]*z[k] > particleRadiusSquared){
						field[i][j][k] = integral(x[i],y[j],z[k],deposits);	
       				}
				currentIteration++;

              			// Display progress bar so that the user has something to look at
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

    //////////////////////////////////////////////////////////////////
    //////////////////////WRITE TO FILE///////////////////////////////
		cout<<"Simulation finished, writing to csv..."<<endl;
		writeFieldToCSV(x,y,z,field);
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


