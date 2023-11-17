#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
#include "functions.h"

double Dt = 
using namespace std;
 	double  bounds;
	double  lambda;
	double  stepSize;
	double  dv;
	int 	nDeposits;	
	int	    nPoints;
	double dl;
	bool onlyTangential = false;

double integral(double x, double y, double z,std::vector<Point> deposits){
	//absorbtionTerm will compute the absorbed ammount of power from the laser
	//ContributionSum will sum up contributions from all deposits
	//Finally, the contributionSum is scaled with volume element dv and divided with constants												
	double laserPower	           = I0 + I0*cos(twoPi*(x)/lambda);	
	double absorbtionTerm          = laserPower*depositArea/(volumePerDeposit);
	double contributionSum 		   = 0.0;
	
	//Since the values scale with the inverse square distance.
    	for (size_t i = 0; i < deposits.size(); i++){

    		double inverse_squareroot_distance = 1.0/sqrt(pow(x-deposits[i].x,2)+
														  pow(y-deposits[i].y,2)+
														  pow(z-deposits[i].z,2));
			contributionSum +=  inverse_squareroot_distance;
		}
    return contributionSum*absorbtionTerm*dv/(4*pi*waterConductivity); 
}
double central_difference(double x1,double y1,double z1,double x2, double y2, double z2, vector<Point> deposits);

int main(int argc, char** argv) {
	   
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    = stold(argv[2])  * pow(10,-6); 
	stepSize  = bounds/(double)(300);		  //Step size, based off of bounds parameter
	nDeposits = stof(argv[1]);				  //number of deposits to initialize
	lambda	  = stold(argv[3])  * pow(10,-9); //Spatial periodicity
    dv	      = stepSize*stepSize*stepSize;  //volume element for integral
	

	int nParticles =1;
	Point centerOfParticle1 = {0.0*particleRadius,0.0,0.0}; 
	Point centerOfParticle2 = {-particleRadius,0.0,0.0}; 
	Particle particle1(centerOfParticle1,particleRadius);
	Particle particle2(centerOfParticle2,particleRadius);
	vector<Particle> particles = {particle1,particle2};

     ///////////GENERATE DEPOSITS//////////////////////////////////////////////////////////////////
	 for(int i = 0; i<nParticles; i++ ){
		particles[i].generateDeposits(nDeposits);
	 }

	//particles[0].rotate(pi/4);
	//particles[1].rotate(-pi/8); 
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
	dl = stepSize*stepSize;
	double thickness = pow(25*stepSize,2);
	double surfaceX = 0.0;
	double surfaceZ = 0.0;
for(int n = 0; n < nParticles;n++){

	
	#pragma omp parallel for
    for (size_t i = 0; i < nPoints; i++){
    	for(size_t j = 0; j<y.size(); j++){
    		for(size_t k = 0; k<nPoints; k++){		
    			//Check if outside particle:
				Point point = {x[i],y[j],z[k]};
				double d = particles[n].getRadialDistance(point);
				if (d > pow(particles[n].radius,2) && d < pow(particles[n].radius,2)+thickness){

					double central_differenceX = central_difference(x[i]-dl,x[i]+dl,y[j],y[j],z[k],z[k],particles[n].deposits);
					
					double central_differenceZ = central_difference(x[i],x[i],y[j],y[j],z[k]-dl,z[k]+dl,particles[n].deposits);
	

					//Project on normal vector:
					double u = x[i]-particles[n].center.x;
					double w = z[k]-particles[n].center.z;

					double perpendicularZ      = (central_differenceX*u+central_differenceZ*w)*w/d;
					double perpendicularX      = (central_differenceX*u+central_differenceZ*w)*u/d;

					//Subtract to get tangential component
					double tangentialX         = (central_differenceX - perpendicularX);
					double tangentialZ         = (central_differenceZ - perpendicularZ);

					if (onlyTangential == true){
						xGrad[i][j][k] 			   += tangentialX*25/1000;		
						zGrad[i][j][k]   		   += tangentialZ*25/1000;			
					}else{
						xGrad[i][j][k] 	     	   += perpendicularX*25/1000;	
						zGrad[i][j][k] 			   += perpendicularZ*25/1000;	
						xGrad[i][j][k] 			   += tangentialX*25/1000;		
						zGrad[i][j][k]   		   += tangentialZ*25/1000;
					}

					surfaceX		  	 += tangentialX*dl*dl;
					surfaceZ			 += tangentialZ*dl*dl;
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
	double total;
	if(onlyTangential == true){
		//Print surface integral (self-propulsion velocity in X and Z components.)
		surfaceX = surfaceX;
		surfaceZ = surfaceZ;
		total = sqrt(pow(surfaceX,2)+pow(surfaceZ,2));
	}
	std::cout<<"Vx = "<< surfaceX<<"\n" << "Vz = "<< surfaceZ<<"\n" << "Vtot = "<<total << "\n";
	//////////////////////////////////////////////////////////////////
	///////////////////COMPUTE ELAPSED TIME///////////////////////////
   	auto endTimer = std::chrono::high_resolution_clock::now();
   	std::chrono::duration<double> duration = endTimer - startTimer;
	double elapsed_seconds = duration.count();
  	std::cout << "Program completed after: " << elapsed_seconds << " seconds" << "\n";
	//////////////////////////////////////////////////////////////////	
	//////////////////END PROGRAM/////////////////////////////////////

	return 0;
}	

double central_difference(double x1,double x2,double y1,double y2, double z1, double z2, vector<Point> deposits){
	double back   		= integral(x1,y1,z1,deposits);
	double forward		= integral(x2,y2,z2,deposits);
	return (forward - back)/(2*dl);
}

