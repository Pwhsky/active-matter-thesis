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

	double f = 1e-9; //This converts to micrometers/s

	//Update positions of deposits and center of particle based on self propulsion
	for(int i = 0; i< this->deposits.size(); i++){
		this->deposits[i].x += (this->selfPropulsion)[0]*f;
		this->deposits[i].z += (this->selfPropulsion)[1]*f;
	}
	this->center.x += (this->selfPropulsion)[0]*f;
	this->center.z += (this->selfPropulsion)[1]*f;

	//Todo: Update particle positions based on external force (from other particles)

}
void Particle::rotate(double angle) {


	//Rotation only works for small angle increments when updating the positions of the deposits
	//during the brownian simulation, the largest possible angle of rotation will be small either way.
	for(int l = 0; l<100;l++){
		double theta =  angle*0.01;

    	for (int i = 0; i < this->deposits.size(); i++) {
       		double distance = getRadialDistance(deposits[i]);
        	this->deposits[i].x = (this->deposits[i].x - this->center.x) * cos(theta) - (this->deposits[i].z - this->center.z) * sin(theta) + this->center.x;
        	this->deposits[i].z = (this->deposits[i].x - this->center.x) * sin(theta) + (this->deposits[i].z - this->center.z) * cos(theta) + this->center.z;
    	}

	}

	double vx = this->selfPropulsion[0];
	double vy = this->selfPropulsion[1];
	double magnitude = sqrt(vx*vx + vy*vy);
	this->selfPropulsion[0] = magnitude*cos(angle);
	this->selfPropulsion[1] = magnitude*sin(angle);
	
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


