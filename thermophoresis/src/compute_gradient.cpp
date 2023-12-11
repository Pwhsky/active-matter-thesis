#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
#include "functions.h"


	constexpr double pi	      		 	  	  = 3.14159265358979323846;
	constexpr double twoPi					  = 2*3.14159265358979323846;
	constexpr double particleRadius   		  = 2    *pow(10,-6);
	constexpr double particleRadiusSquared    = particleRadius*particleRadius;
	constexpr double depositRadius	 	      = 30   *pow(10,-9);	
	constexpr double volumePerDeposit	      = 4*pi *pow((depositRadius),3)/3; 
	constexpr double depositArea	 	      = 2*pi *pow(depositRadius,2); 
	constexpr double intensity		  		  = 100  *pow(10,-3);  //Watts
	constexpr double areaOfIllumination 	  = 40   *pow(10,-6);  //Meters  How much area the laser is distributed on.
	constexpr double I0		 				  = 2*intensity/(pow(areaOfIllumination*2,2)); 
	constexpr double waterConductivity	 	  = 0.606;


using namespace std;
 	double  bounds;
	double  lambda;
	double  stepSize;
	double  dv;
	int 	nDeposits;	
	int	    nPoints;
	double dl;
	double thickness; 
	bool onlyTangential = false;

double thermoDiffusion = 2.8107e-6; 
void getSelfPropulsion(Particle particle,vector<double> x, 
						vector<double> y,vector<double> z);


double integral(double x, double y, double z,std::vector<Point> deposits);
double central_difference(double x1,double y1,
						  double z1,double x2, 
						  double y2, double z2,
						vector<Point> deposits);

int main(int argc, char** argv) {
	   
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    = stold(argv[2])  * pow(10,-6); 
	stepSize  = bounds/(double)(300);		  //Step size, based off of bounds parameter
	nDeposits = stof(argv[1]);				  //number of deposits to initialize
	lambda	  = stold(argv[3])  * pow(10,-9); //Spatial periodicity
    dv	      = stepSize*stepSize*stepSize;  //volume element for integral
	

	

	vector<Particle> particles = initializeParticles();

	int nParticles = particles.size();
	for(int i = 0; i<nParticles; i++ )
		particles[i].generateDeposits(nDeposits);

	//particles[0].rotate(pi/4);
	//particles[1].rotate(-pi/8); 
	//////////////////////////////////////////////////////////////////////////////////////////////
	///////////INITIALIZE LINSPACE VECTORS////////////////////////////////////////////////////////
	vector<double> linspace = arange(-bounds,bounds,stepSize);
	vector<double> y = {0.0};
    vector<double> z = linspace;
	vector<double> x = linspace;
	nPoints = x.size();
   
 	cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
	vector<vector<vector<double>>> xGrad(nPoints, vector<vector<double>>(nPoints, vector<double>(nPoints)));
    vector<vector<vector<double>>> zGrad(nPoints, vector<vector<double>>(nPoints, vector<double>(nPoints)));
      
    const int totalIterations = nPoints*nPoints*y.size();
	size_t currentIteration = 0;


	dl = stepSize*stepSize;
	thickness = pow(25*stepSize,2);

	for(int n = 0; n < nParticles;n++){

		getSelfPropulsion(particles[n],x,y,z);
		double Qx = particles[n].center.x;
		double Qy = particles[n].center.y;
		double Qz = particles[n].center.z;


		int counter = 0;
		#pragma omp parallel for
		for (size_t i = 0; i < nPoints; i++){

			for(size_t j = 0; j<y.size(); j++){
				for(size_t k = 0; k<nPoints; k++){		
					//Check if outside particle:
					Point point = {x[i],y[j],z[k]};
					double d = particles[n].getRadialDistance(point);
					if (d > pow(particles[n].radius,2)){

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

void getSelfPropulsion(Particle particle, vector<double> x,vector<double> y,vector<double> z){
	double Qx = particle.center.x;
	double Qy = particle.center.y;
	double Qz = particle.center.z;
	int counter = 1;

	#pragma omp parallel for
		for (size_t i = 0; i < nPoints; i++){
			for(size_t j = 0; j<y.size(); j++){
				for(size_t k = 0; k<nPoints; k++){		
					
					Point point = {x[i],y[j],z[k]};
					double d = particle.getRadialDistance(point);
					//Compute only the points near the surface
					if (d > pow(particle.radius,2) && d < pow(particle.radius,2)+thickness){

						double central_differenceX = central_difference(x[i]-dl,x[i]+dl,y[j],y[j],z[k],z[k],particle.deposits);
						double central_differenceZ = central_difference(x[i],x[i],y[j],y[j],z[k]-dl,z[k]+dl,particle.deposits);

						//Project on normal vector:
						double u = x[i]-particle.center.x;
						double w = z[k]-particle.center.z;

						double perpendicularZ      = (central_differenceX*u+central_differenceZ*w)*w/d;
						double perpendicularX      = (central_differenceX*u+central_differenceZ*w)*u/d;

						//Subtract to get tangential component
						double tangentialX         = (central_differenceX - perpendicularX);
						double tangentialZ         = (central_differenceZ - perpendicularZ);
		
		
						//Surface integral to compute self-propulsion

						double theta = atan2(sqrt((Qx-x[i])*(Qx-x[i]) + (Qz-z[k])*(Qz-z[k])), sqrt((Qy-y[j]) *(Qy-y[j]) ));
						
						double phi   = atan((Qz-z[k])/(Qx-x[i]));
						particle.selfPropulsion[0] += tangentialX*sin(theta)*cos(phi);
						particle.selfPropulsion[1] += tangentialZ*sin(theta)*cos(phi);
						counter++;
					}
					
				}
			}
		}
		//Scale with diffusion coefficient
	for(int i = 0; i<1;i++){
		particle.selfPropulsion[i] *= pi*pi*thermoDiffusion/(double)counter;
		//cout<<"\n"<<particle.selfPropulsion[i]<<"\n";
	}
		
}
