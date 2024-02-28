/*
Alex Lech 2023	

This file contains the functions source code.

In functions.h the particle class is declared.

*/


#include "particle.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <random>

using namespace std;
std::random_device rd;
std::mt19937 gen(rd());

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
	constexpr long double dt = 0.1; 
	
	

void Particle::generateDeposits(int nDeposits) {
	//Create random number generators, with certain intervals.
	//to create a janus-particle distribution, costheta parameter should be between 0 and 1 (corresponding to z axis).
	//If one wants a completely covered particle, set costheta to (-1,1).
	//To adjust the areas of initialization, play around with the phi and costheta parameters.
	//u is commonly set (0.9,1) so the deposits are near the surface.

    uniform_real_distribution<double> phi(0.0,twoPi); 
    uniform_real_distribution<double> costheta(-0.1,1);
    uniform_real_distribution<double> u(0.8,1);

	//Initiate deposits
	for(int i = 0; i<nDeposits; i++){

		double theta = acos(costheta(gen));
		double r 	 = (this->radius)*u(gen);

		//Convert to cartesian:
    	double x = r*sin(theta) * cos(phi(gen)) + this->center.x; 
    	double y = r*sin(theta) * sin(phi(gen)) + this->center.y;
    	double z = r*cos(theta)					+ this->center.z;
   		
		//Add to deposits list
    	(this->deposits).emplace_back(Point{x,y,z});
    	
    }

	//Some interresting configurations:
	//phi(0.0,pi/2)
	//costheta(0.5,0.7)
	//u(0.9,1)

	//phi(0.0,pi/4)
	//costheta(0.5,0.7)
	//u(0.9,1)
}
double Particle::getRadialDistance(Point r){
		double norm = pow(this->center.x-r.x,2) + pow(this->center.y-r.y,2) + pow(this->center.z-r.z,2);
	   return norm;
}

void Particle::updatePosition(){

	//Update positions of deposits and center of particle based on self propulsion
	for(int i = 0; i< this->deposits.size(); i++){

		this->deposits[i].x += (this->velocity)[0]*dt;
		this->deposits[i].y += (this->velocity)[1]*dt;
		this->deposits[i].z += (this->velocity)[2]*dt;
	}
	this->center.x += (this->velocity)[0]*dt;
	this->center.y += (this->velocity)[1]*dt;
	this->center.z += (this->velocity)[2]*dt;

}

//      This function extensively uses the following formulae for rotation: 
//          Numerical Simulations of Active Brownian Particles by A. Callegari & G. Volpe
//          Section 7.4.1

void Particle::rotation_transform() {
    double* w = this->selfRotation;
    double theta = dt* sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);

    if (theta != 0) {

        //Generate theta_x matrix
        std::vector<std::vector<double>> theta_x = {
            { 0, -w[2], w[1] },
            { w[2], 0, -w[0] },
            { -w[1], w[0], 0 }
        };

        std::vector<std::vector<double>> theta_x_squared = mat_mat_mul(theta_x, theta_x);
        std::vector<std::vector<double>> R = theta_x;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                R[i][j] = (sin(theta) / theta) * theta_x[i][j] + ((1.0 - cos(theta)) / (theta * theta)) * theta_x_squared[i][j];

                if (i == j){
                    R[i][j] +=1;
                }
            }
			
        }

        // Rotate the deposits, using the particle center as reference

        std::vector<double> temp(3);
        std::vector<double> x(3);
   
        for (auto &p : this->deposits) {

            x = { p.x-this->center.x, p.y-this->center.y, p.z-this->center.z };
            temp = mat_vec_mul(R,x);

            p.x = temp[0];
            p.y = temp[1];
            p.z = temp[2];
        }
    
        
        // Rotate the particle itself (center)
        x = { this->center.x, this->center.y, this->center.z };
        std::vector<double> v(this->velocity, this->velocity +3);
        temp = mat_vec_mul(R,v);
        for(int i  = 0; i<3; i++){
            this->velocity[i] = temp[i];
        }

    }
}


double integral(double _x,double _y,double _z,std::vector<Point> deposits,double lambda,double dv){
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

double central_difference(double x1,double x2,double y1,double y2, double z1, double z2, vector<Point> deposits,double dl,double lambda, double dv){
	double back   		= integral(x1,y1,z1,deposits,lambda,dv);
	double forward		= integral(x2,y2,z2,deposits,lambda,dv);
	return (forward - back)/(2*dl);
}




void Particle::getKinematics(std::vector<double> x,std::vector<double> y, std::vector<double> z,
				double thickness,double dl,std::vector<Point> globalDeposits, double lambda, double dv){

	//This will compute the tangential component in a thin layer around the particle
	//And then do a surface integral to get self propulsion in X and Z direction.


	int counter = 1;
	vector<double> omega = {0,0,0};

	double gradientX;
	double gradientY;
	double gradientZ;
	double Qx = center.x;
	double Qy = center.y;
	double Qz = center.z;

	#pragma omp parallel for

		for (auto i:x){
			for(auto j:y){
				for(auto k:z){		
					
					Point point 	= { i-center.x,
										j-center.y,
										k-center.z};

					//double norm = get_norm({point.x, point.y, point.z});
					//if (norm > 6*particleRadius){continue;}

					double d = getRadialDistance(point);


					//Compute only the points near the surface
					if (d > pow(radius,2) && d < pow(radius,2)+thickness){
							
							double u				   = point.x;
							double v 				   = point.y;
							double w 				   = point.z;

							vector<double> r		   = {u,v,w};	


							vector<double> gradient = {central_difference(i-dl,i+dl,j   ,j   ,k   ,k   ,globalDeposits,dl,lambda,dv),
													   central_difference(i   ,i   ,j-dl,j+dl,k   ,k   ,globalDeposits,dl,lambda,dv),
													   central_difference(i   ,i   ,j   ,j   ,k-dl,k+dl,globalDeposits,dl,lambda,dv)};

	
							vector<double> radial     = {0,0,0};
							vector<double> tangential = {0,0,0};
							vector<double> vel 		  = {0,0,0};
							double duvwr = gradient[0] * u + gradient[1] * v + gradient[2] * w;
							double theta = atan2(sqrt((Qx-i)*(Qx-i) + (Qz-k)*(Qz-k)), sqrt((Qy-j) *(Qy-j) ));
							double phi   = atan((Qz-k)/(Qx-i));
							double sincos = sin(theta)*cos(phi);

							//Populate the cartesian vectors with their respective components:
							for(int l = 0; l<3; l++){
								radial[l]    	   = duvwr * r[l] / d;
								tangential[l]	   = gradient[l] - radial[l];
								vel[l]     		   = tangential[l] * sincos;
								velocity[l] 	   = vel[l];
							}

							vector<double> rxV = cross_product(r,vel);

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
		velocity[i]    *= thermoDiffusion/(double)counter;
		selfRotation[i] = omega[i]/((double)counter);
		//cout<<output.selfRotation[i]<<"\n";
		//cout<<output.velocity[i]<<"\n";
	}

}


void Particle::writeDepositToCSV() {
    static bool isFirstRun = true;

    std::ofstream outputFile;

    if (isFirstRun) {
        // If it's the first run, create a new file with the header
        outputFile.open("deposits.csv");
        outputFile << "x,y,z" << "\n";
        isFirstRun = false;
    } else {
        // If it's not the first run, open the file in append mode
        outputFile.open("deposits.csv", std::ios::app);
    }

    // Write data to the file
    for (size_t i = 0; i < size(deposits); i++) {
        outputFile << (this->deposits)[i].x << "," << (this->deposits)[i].y  << "," << (this->deposits)[i].z  << "\n";
    }

    outputFile.close();
}





//Old rotation function
/*
void Particle::rotate(double angle) {
	//Rotation only works for small angle increments when updating the positions of the deposits
	//during the brownian simulation, the largest possible angle of rotation will be small either way.
	for(int l = 0; l<100;l++){
		double theta =   angle*0.01;

    	for (int i = 0; i < this->deposits.size(); i++) {
       		double distance = getRadialDistance(deposits[i]);
        	this->deposits[i].x = (this->deposits[i].x - this->center.x) * cos(theta) - (this->deposits[i].z - this->center.z) * sin(theta) + this->center.x;
        	this->deposits[i].z = (this->deposits[i].x - this->center.x) * sin(theta) + (this->deposits[i].z - this->center.z) * cos(theta) + this->center.z;
    	}

	}

	double vx = this->velocity[0];
	double vy = this->velocity[1];
	double vz = this->velocity[2];

	double magnitude = sqrt(vx*vx + vy*vy + vz*vz);

	this->velocity[0] = magnitude*cos(angle);
	this->velocity[1] = magnitude*sin(angle);
	
}
*/



