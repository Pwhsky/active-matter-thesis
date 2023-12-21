/*
Alex Lech 2023	

This source code contains the computation for the temperature gradient.
compute_temperature.cpp contains the computation for the temperature increase.
*/

#include <chrono>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include "particle.h"
#include <random>
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
	constexpr double thermoDiffusion 		  = 2.8107e-6; 


using namespace std;
 	double  bounds;
	double  lambda;
	double  stepSize;
	double  dv;
	int 	nDeposits;	
	int	    nPoints;
	double  dl;
	double  thickness; 
	bool    onlyTangential = false;
	vector<double> z,x,y;

Particle getSelfPropulsion(Particle particle);
double integral(double _x, double _y,double _z,std::vector<Point> deposits);

double central_difference(double x_back,double x_forward,
						  double y_back,double y_forward, 
						  double z_back, double z_forward,
						  vector<Point> deposits);

int main(int argc, char** argv) {
	   
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    = stold(argv[2])  * pow(10,-6); 
	stepSize  = bounds/(double)(300);		  //Step size, based off of bounds parameter
	nDeposits = stof(argv[1]);				  //number of deposits to initialize
	lambda	  = stold(argv[3])  * pow(10,-9); //Spatial periodicity
    dv	      = stepSize*stepSize*stepSize;  //volume element for integral

	
	std::ofstream writePositions("positions.csv");
	writePositions << "x,y,z"<<"\n";

	vector<Particle> particles = initializeParticles();
	int nParticles = particles.size();
	for(int i = 0; i<nParticles; i++ ) particles[i].generateDeposits(nDeposits);

	//Generate linspace vectors:
	vector<double> linspace = arange(-bounds,bounds,stepSize);
	z = linspace;
	x = linspace;
	y = linspace;

	nPoints = x.size();
   
 	cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
	vector<vector<vector<double>>> xGrad(nPoints, vector<vector<double>>(nPoints, vector<double>(nPoints)));
    vector<vector<vector<double>>> zGrad(nPoints, vector<vector<double>>(nPoints, vector<double>(nPoints)));
      
    const int totalIterations = nPoints*nPoints*y.size();
	size_t currentIteration;

	dl = stepSize*stepSize;
	thickness = pow(10*stepSize,2);

	
	for(int n = 0; n<nParticles;n++){
		particles[0] = getSelfPropulsion(particles[0]);
	}
	std::random_device rd;
	std::mt19937 gen(rd());
	
	uniform_real_distribution<double> phi(-pi,pi); 
	for(int time = 0; time < 100; time ++){
		currentIteration = 0;

		for(int n = 0; n < nParticles;n++){
			
			
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
							/*
							
							
							double gradientX = central_difference(x[i]-dl,x[i]+dl,y[j],   y[j],   z[k],   z[k],   particles[n].deposits);
							double gradientY = central_difference(x[i],   x[i],   y[j]-dl,y[j]+dl,z[k],   z[k],   particles[n].deposits);
							double gradientZ = central_difference(x[i],   x[i],   y[j],   y[j],   z[k]-dl,z[k]+dl,particles[n].deposits);

							//Project on normal vector:
							double u 				   = x[i]-Qx;
							double v				   = y[j]-Qy;
							double w 				   = z[k]-Qz;
							

							double perpendicularX      = (gradientX*u + gradientY*v + gradientZ*w)*u/d;
							double perpendicularY      = (gradientX*u + gradientY*v + gradientZ*w)*v/d;
							double perpendicularZ      = (gradientX*u + gradientY*v + gradientZ*w)*w/d;

							//Subtract to get tangential component
							double tangentialX         = (gradientX - perpendicularX);
							double tangentialY		   = (gradientY - perpendicularY);
							double tangentialZ         = (gradientZ - perpendicularZ);
							
							
						
							
							if (onlyTangential == true){
								xGrad[i][j][k] 			   += tangentialX*25/1000;		
								zGrad[i][j][k]   		   += tangentialZ*25/1000;			
							}else{
								xGrad[i][j][k] 	     	   += perpendicularX*25/1000;	
								xGrad[i][j][k] 			   += tangentialX*25/1000;	

								zGrad[i][j][k] 			   += perpendicularZ*25/1000;	
								zGrad[i][j][k]   		   += tangentialZ*25/1000;
							}

							*/
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
			//particles[n].rotate(phi(gen));
			particles[n].updatePosition();
			
			
			writePositions << particles[n].center.x << "," << particles[n].center.y << "," << particles[n].center.z << "\n";
			//particles[n].writeDepositToCSV();
		}
		cout<<"Finished step "<<time<<"\n";
	}
	particles[0].writeDepositToCSV();

    //////////////////////////////////////////////////////////////////
    //////////////////////WRITE TO FILE///////////////////////////////
	cout<<"Simulation finished, writing to csv..."<<"\n";
	//writeGradToCSV(x,y,z,xGrad,zGrad);
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

inline double central_difference(double x_back,double x_forward,double y_back,double y_forward, double z_back, double z_forward, vector<Point> deposits){
	double back   		= integral(x_back,y_back,z_back,deposits);
	double forward		= integral(x_forward,y_forward,z_forward,deposits);
	return (forward - back)/(2*dl);
}

inline double integral(double _x,double _y,double _z,std::vector<Point> deposits){
	//absorbtionTerm will compute the absorbed ammount of power from the laser
	//ContributionSum will sum up contributions from all deposits
	//Finally, the contributionSum is scaled with volume element dv and divided with constants												
	double laserPower	           = I0 + I0*cos(twoPi*(_x)/lambda);	
	double absorbtionTerm          = laserPower*depositArea/(volumePerDeposit);
	double contributionSum 		   = 0.0;
	
	//Since the values scale with the inverse square distance.
    	for (size_t i = 0; i < deposits.size(); i++){

    		double inverse_squareroot_distance = 1.0/sqrt(pow(_x-deposits[i].x,2)+
														  pow(_y-deposits[i].y,2)+
														  pow(_z-deposits[i].z,2));
			contributionSum +=  inverse_squareroot_distance;
		}
    return contributionSum*absorbtionTerm*dv/(4*pi*waterConductivity); 
}

Particle getSelfPropulsion(Particle particle){
	//This will compute the tangential component in a thin layer around the particle
	//And then do a surface integral to get self propulsion in X and Z direction.

	double Qx = particle.center.x;
	double Qy = particle.center.y;
	double Qz = particle.center.z;
	int counter = 1;

		#pragma omp parallel for
		for (auto i:x){
			for(auto j:y){
				for(auto k:z){		
					
					Point point = {i,j,k};
					double d = particle.getRadialDistance(point);
					//Compute only the points near the surface
					if (d > pow(particle.radius,2) && d < pow(particle.radius,2)+thickness){
						
						double gradientX = central_difference(i-dl,i+dl,j,j,k,k,particle.deposits);
						double gradientY = central_difference(i,i,j-dl,j+dl,k,k,particle.deposits);
						double gradientZ = central_difference(i,i,j,j,k-dl,k+dl,particle.deposits);

						//Project on normal vector:
						double u				   = i-particle.center.x;
						double v 				   = j-particle.center.y;
						double w 				   = k-particle.center.z;

						double perpendicularZ      = (gradientX*u + gradientY*v + gradientZ*w)*w/d;
						double perpendicularY      = (gradientX*u + gradientY*v + gradientZ*w)*v/d;
						double perpendicularX      = (gradientX*u + gradientY*v + gradientZ*w)*u/d;

						//Subtract to get tangential component
						double tangentialX         = (gradientX - perpendicularX);
						double tangentialY         = (gradientY - perpendicularY);
						double tangentialZ         = (gradientZ - perpendicularZ);

						double theta = atan2(sqrt((Qx-i)*(Qx-i) + (Qz-k)*(Qz-k)), sqrt((Qy-j) *(Qy-j) ));
						double phi   = atan((Qz-k)/(Qx-i));
						particle.selfPropulsion[0] += tangentialX*sin(theta)*cos(phi);
						particle.selfPropulsion[1] += tangentialY*sin(theta)*cos(phi);
						particle.selfPropulsion[2] += tangentialZ*sin(theta)*cos(phi);
						counter++;
					}
					
				}
			}
		}

		
		//Scale with diffusion coefficient
	for(int i = 0; i<3;i++){
		particle.selfPropulsion[i] *= thermoDiffusion/(double)counter;
	
		//cout<<"\n"<<particle.selfPropulsion[i]<<"\n";
	}

	return particle;		
}
