#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
#include "functions.h"

using namespace std;

	long double stepSize;
	long double dv;


inline  long double integral(long double x, long double y, long double z,vector<Point> deposits){
	
	long double laserTerm 	    = I0 + I0*cos(twoPi*x/lambda);
	long double qTerm 	    = laserTerm*depositArea/(volumePerDeposit);
	long double contributionSum    = 0.0;
	long double q                  = 0.0;


    	for (size_t i = 0; i < deposits.size(); i++){
    		long double inv_sqrt_distance1 = 1.0/sqrt((x-deposits[i].x)*(x-deposits[i].x) + 
    					      (y-deposits[i].y)*(y-deposits[i].y) +
    					      (z-deposits[i].z)*(z-deposits[i].z));
  		
		contributionSum +=  inv_sqrt_distance1;
	}
    	return contributionSum*dv*qTerm/(4*pi*waterConductivity);
}




int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();

	stepSize = bounds/(stof(argv[1]));
        dv       = stepSize*stepSize*stepSize; //volume element for integral


        ///////////GENERATE DEPOSITS//////////////////////////////////////////////////////////////////
	vector<Point> deposits;
	generateDeposits(deposits);
	//////////////////////////////////////////////////////////////////////////////////////////////
	///////////INITIALIZE LINSPACE VECTORS////////////////////////////////////////////////////////
	 vector< long double> spaceVector;
	 for (long double coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
          	 spaceVector.push_back(coordinate);
      	 }
     	 const vector<long double> z = spaceVector;
	 const vector<long double> x = spaceVector;
	 vector<long double>       y = {0.0};
	 vector<vector<vector<long double>>> field(x.size(), vector<vector<long double>>(y.size(), vector<long double>(z.size())));      
 	 cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
 	//////////////////////////////////////////////////////////////////////////////////////////////
 	////////////////////////  INTEGRAL ///////////////////////////////////////////////////////////
	const int totalIterations = x.size() * y.size() * z.size();
   	int currentIteration = 0;
 	#pragma omp parallel for
    	for (int i = 0; i<x.size(); i++){
    		for(int j = 0; j<y.size(); j++){
    			for(int k = 0; k<z.size(); k++){

				field[i][j][k] = integral(x[i],y[j],z[k],deposits);
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


