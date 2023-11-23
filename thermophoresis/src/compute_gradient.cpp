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
	double  dl;

	//When only computing tangential flow and self-propulsion:
	bool onlyTangential = true;

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

	double velocity = 0.0;
	double thermoDiffusion = 2.8107e-6; //meters^2 / kelvin*second
	Particle particle1(centerOfParticle1,particleRadius,velocity);
	Particle particle2(centerOfParticle2,particleRadius,velocity);
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

	double surfaceX = 0.0;
	double surfaceZ = 0.0;
	double counter = 0;
	double thickness = pow(15*stepSize,2);

	//X and Z components of the gradient
	double perpendicular[2];
	double tangential[2];
	double central_differences[2];
	double surface_integral[2];

	for(int n = 0; n < nParticles;n++){
		#pragma omp parallel for
		for (size_t i = 0; i < nPoints; i++){
			for(size_t j = 0; j<y.size(); j++){
				for(size_t k = 0; k<nPoints; k++){	
						
					//Check if outside particle:
					Point point = {x[i],y[j],z[k]};
					double d = particles[n].getRadialDistance(point);
					if (d > pow(particles[n].radius,2) && d < pow(particles[n].radius,2)+thickness){

						central_differences[0] = central_difference(x[i]-dl,x[i]+dl,y[j],y[j],z[k],   z[k],   particles[n].deposits);
						central_differences[1] = central_difference(x[i],   x[i],   y[j],y[j],z[k]-dl,z[k]+dl,particles[n].deposits);
		

						//Project on normal vector:
						double u = x[i]-particles[n].center.x; //u
						double w = z[k]-particles[n].center.z; //w

						double gradientSum = (central_differences[0]*u + central_differences[1]*w)/d;

						
						perpendicular[0] =  u*gradientSum;
						perpendicular[1] =  w*gradientSum;

						//Subtract to get tangential component
						tangential[0]	 = central_differences[0] - perpendicular[0];
						tangential[1]	 = central_differences[1] - perpendicular[1];

						xGrad[i][j][k] 			   += tangential[0]*25/1000.0;		
						zGrad[i][j][k]   		   += tangential[1]*25/1000.0;

						if (onlyTangential == false){
							xGrad[i][j][k] 	     	   += perpendicular[0]*25/1000.0;	
							zGrad[i][j][k] 			   += perpendicular[1]*25/1000.0;				
						}

						//Surface integral to compute self-propulsion

						double theta = atan2(sqrt(x[i]*x[i] + z[k]*z[k]), y[j]);
						double phi   = atan(z[k]/x[i]);
						particles[n].selfPropulsion[0] += tangential[0]*sin(theta)*cos(phi);
						particles[n].selfPropulsion[1] += tangential[1]*sin(theta)*cos(phi);

						counter++;
					}
					currentIteration++;


					//Progress bar
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
		
		//scale the values with the thermodiffusion coefficient and the number of points computed.
		for(int t = 0; t<2; t++)
			particles[n].selfPropulsion[t] *= pi*pi*thermoDiffusion/(double)counter;
		
		
	}


	cout<<"Simulation finished, writing to csv..."<<endl;
	writeGradToCSV(x,y,z,xGrad,zGrad);



	std::cout<<"Vx = "<< particles[0].selfPropulsion[0]<<"\n" << "Vz = "<< particles[0].selfPropulsion[1]<<"\n";
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

