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
	int     number_of_steps,nDeposits, nPoints;
 	double  bounds, lambda, stepSize, dv, dl, thickness; 
	bool    onlyTangential = false;
	vector<double> z,x,y;
	vector<Point> globalDeposits;


inline Particle getKinematics(Particle particle);
inline double integral(double _x, double _y,double _z,std::vector<Point> deposits);

inline double central_difference(double x_back,double x_forward,
						  		double y_back,double y_forward, 
						  		double z_back, double z_forward,
						  		vector<Point> deposits);


int main(int argc, char** argv) {
	   
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    		= stold(argv[2])  * pow(10,-6); 
	stepSize  		= bounds/(double)(300);		  //Step size, based off of bounds parameter Set to 100 for now
	nDeposits 		= stof(argv[1]);				  //number of deposits to initialize
	lambda	 	    = stold(argv[3])  * pow(10,-9); //Spatial periodicity
	number_of_steps = (int)stof(argv[4]);
    dv	      		= stepSize*stepSize*stepSize;  //volume element for integral

	
	std::ofstream writeInitial("initpositions.csv");
	std::ofstream writeFinal("finalpositions.csv");

	writeInitial << "x,y,z"<<"\n";
	writeFinal   << "x,y,z"<<"\n";

	vector<Particle> particles = initializeParticles();
	int nParticles = particles.size();

	for(int i = 0; i<nParticles; i++ ){
		particles[i].generateDeposits(nDeposits);
		for(int j = 0; j<nDeposits;j++){
			globalDeposits.push_back(particles[i].deposits[j]);
		}
	} 


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

	dl = stepSize*stepSize;
	thickness = pow(10*stepSize,2);

	
	//Write initial positions:
	for(auto p:particles) writeInitial << p.center.x << "," << p.center.y << "," << p.center.z << "\n";
	
	
	for(int time = 0; time < number_of_steps; time ++){ //Loop over timestep

		for(auto &particle:particles){
			particle = getKinematics(particle);
			particle.updatePosition();
			particle.rotation_transform();
		}
			
		cout<<"Finished step "<<time<<"\n";
	}

	//Write final position
	for(auto p:particles) writeFinal << p.center.x << "," << p.center.y << "," << p.center.z << "\n";
	

	particles[0].writeDepositToCSV();
	particles[1].writeDepositToCSV();

	cout<<"Simulation finished, writing to csv..."<<"\n";
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

Particle getKinematics(Particle particle){
	//This will compute the tangential component in a thin layer around the particle
	//And then do a surface integral to get self propulsion in X and Z direction.

	Particle output = particle;

	int counter = 1;
	vector<double> omega = {0,0,0};

	double gradientX;
	double gradientY;
	double gradientZ;
	double Qx = particle.center.x;
	double Qy = particle.center.y;
	double Qz = particle.center.z;

	#pragma omp parallel for

		for (auto i:x){
			for(auto j:y){
				for(auto k:z){		
					
					Point point 	= { i-particle.center.x,
										j-particle.center.y,
										k-particle.center.z};

					//double norm = get_norm({point.x, point.y, point.z});
					//if (norm > 6*particleRadius){continue;}

					double d = particle.getRadialDistance(point);


					//Compute only the points near the surface
					if (d > pow(particle.radius,2) && d < pow(particle.radius,2)+thickness){
							
							double u				   = point.x;
							double v 				   = point.y;
							double w 				   = point.z;

							vector<double> r		   = {u,v,w};	


							vector<double> gradient = {central_difference(i-dl,i+dl,j   ,j   ,k   ,k   ,globalDeposits),
													   central_difference(i   ,i   ,j-dl,j+dl,k   ,k   ,globalDeposits),
													   central_difference(i   ,i   ,j   ,j   ,k-dl,k+dl,globalDeposits)};

	
							vector<double> radial     = {0,0,0};
							vector<double> tangential = {0,0,0};
							vector<double> velocity   = {0,0,0};

							double duvwr = gradient[0] * u + gradient[1] * v + gradient[2] * w;
							double theta = atan2(sqrt((Qx-i)*(Qx-i) + (Qz-k)*(Qz-k)), sqrt((Qy-j) *(Qy-j) ));
							double phi   = atan((Qz-k)/(Qx-i));
							double sincos = sin(theta)*cos(phi);

							//Populate the cartesian vectors with their respective components:
							for(int l = 0; l<3; l++){
								radial[l]    	   = duvwr * r[l] / d;
								tangential[l]	   = gradient[l] - radial[l];
								velocity[l]        = tangential[l] * sincos;
								output.velocity[l] = velocity[l];
							}

							vector<double> rxV = cross_product(r,velocity);

							for(int l = 0; l<3; l++){
								omega[l] -= rxV[l];
							}

							counter++;

					}
					
				}
			}
		}
		//scale with number of points
	for(int i = 0; i<3;i++){
		output.velocity[i]    *= thermoDiffusion/(double)counter;
		output.selfRotation[i] = omega[i]/((double)counter);
		cout<<output.selfRotation[i]<<"\n";
		//cout<<output.velocity[i]<<"\n";
	}

	return output;		
}


/*

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
						
						double gradientX = central_difference(i-dl,i+dl,j   ,j   ,k   ,k   ,globalDeposits);
						double gradientY = central_difference(i   ,i   ,j-dl,j+dl,k   ,k,globalDeposits);
						double gradientZ = central_difference(i   ,i   ,j   ,j   ,k-dl,k+dl,globalDeposits);

						//Project on normal vector:
						double u				   = i-particle.center.x;
						double v 				   = j-particle.center.y;
						double w 				   = k-particle.center.z;

						double perpendicularX      = (gradientX*u + gradientY*v + gradientZ*w)*u/d;
						double perpendicularY      = (gradientX*u + gradientY*v + gradientZ*w)*v/d;
						double perpendicularZ      = (gradientX*u + gradientY*v + gradientZ*w)*w/d;

						//Subtract to get tangential component
						double tangentialX         = (gradientX - perpendicularX);
						double tangentialY         = (gradientY - perpendicularY);
						double tangentialZ         = (gradientZ - perpendicularZ);

						double theta = atan2(sqrt((Qx-i)*(Qx-i) + (Qz-k)*(Qz-k)), sqrt((Qy-j) *(Qy-j) ));
						double phi   = atan((Qz-k)/(Qx-i));

						particle.velocity[0] += tangentialX*sin(theta)*cos(phi);
						particle.velocity[1] += tangentialY*sin(theta)*cos(phi);
						particle.velocity[2] += tangentialZ*sin(theta)*cos(phi);
						counter++;
					}
					
				}
			}
		}
		
		//Rescale:
	for(int i = 0; i<3;i++){
		particle.velocity[i] *= thermoDiffusion/(double)counter;
	}

	return particle;		
}

	
			
	 		int counter = 0;
			double Qx = particles[n].center.x;
			double Qy = particles[n].center.y;
			double Qz = particles[n].center.z;

			#pragma omp parallel for
			for (size_t i = 0; i < nPoints; i++){
				for(size_t j = 0; j<y.size(); j++){
					for(size_t k = 0; k<nPoints; k++){

						//Check if outside particle:

						Point point = {x[i],y[j],z[k]};
						double d = particles[n].getRadialDistance(point);
						if (d > pow(particles[n].radius,2)){
							
						
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
						}
						

						// Calculate progress percentage so that the user has something to look at
						currentIteration++;
						if(currentIteration % 10000 == 0) {
							float progress = round(static_cast<float>(currentIteration) / (nParticles*totalIterations) * 100.0);
							// Print progress bar
							#pragma omp critical
							{
									cout << "Progress: "<< progress << "% \r";
									cout.flush();
							}
						}
					}
				}
			}
			*/