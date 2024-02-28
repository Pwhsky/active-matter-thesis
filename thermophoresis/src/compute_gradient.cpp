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

	vector<double> z;
	vector<double> x;
	vector<double> y;


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
	y = {0.0};

	nPoints = x.size();
   
 	cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
	vector<vector<vector<double>>> xGrad(nPoints, vector<vector<double>>(nPoints, vector<double>(nPoints)));
    vector<vector<vector<double>>> zGrad(nPoints, vector<vector<double>>(nPoints, vector<double>(nPoints)));
      
    const int totalIterations = nPoints*nPoints*y.size();
	size_t currentIteration = 0;


	dl = stepSize*stepSize;
	thickness = pow(10*stepSize,2);


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

							double gradientX = central_difference(x[i]-dl,x[i]+dl,y[j],y[j],z[k],   z[k],   particles[n].deposits,dl,lambda,dv);
							double gradientZ = central_difference(x[i],   x[i],   y[j],y[j],z[k]-dl,z[k]+dl,particles[n].deposits,dl,lambda,dv);
			

							//Project on normal vector:
							double u 				   = x[i]-Qx;
							double w 				   = z[k]-Qz;
							
							double perpendicularX      = (gradientX*u+gradientZ*w)*u/d;
							double perpendicularZ      = (gradientX*u+gradientZ*w)*w/d;

							//Subtract to get tangential component
							double tangentialX         = (gradientX - perpendicularX);
							double tangentialZ         = (gradientZ - perpendicularZ);

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
						currentIteration++

					}
				}
			}
			
			writePositions << particles[n].center.x << "," << particles[n].center.y << "," << particles[n].center.z << "\n";
		}
		cout<<"Finished step "<<time<<"\n";
	
	particles[0].writeDepositToCSV();

    //////////////////////////////////////////////////////////////////
    //////////////////////WRITE TO FILE///////////////////////////////
	cout<<"Simulation finished, writing to csv..."<<"\n";
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



