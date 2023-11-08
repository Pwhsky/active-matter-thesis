#include "functions.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <random>

using namespace std;
std::random_device rd;
std::mt19937 gen(rd());

	

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

//Work in progress
void Particle::rotate(double theta){
	for(int i = 0;i< (this->deposits).size();i++){
		double distance = getRadialDistance(deposits[i]);
	
		 (this->deposits)[i].x = this->center.x -cos(theta) * (this->deposits)[i].x - sin(theta) * (this->deposits)[i].z;
        (this->deposits)[i].z = this->center.z - sin(theta) * (this->deposits)[i].x + cos(theta) * (this->deposits)[i].z;
	}
}

void Particle::writeDepositToCSV() {
	std::ofstream outputFile("deposits.csv");
    outputFile << "x,y,z" << "\n";
	for (size_t i = 0; i < size(deposits); i++) {
	        outputFile << (this->deposits)[i].x << "," << (this->deposits)[i].y  << "," << (this->deposits)[i].z  << "\n";
	}

    outputFile.close();
}

void writeGradToCSV(const std::vector<double>& x, 
		     const std::vector<double>& y, 
		     const std::vector<double>& z, 
		     std::vector<std::vector<std::vector<double>>>& xGrad,
			 std::vector<std::vector<std::vector<double>>>& yGrad){

	std::ofstream outputFile("gradient.csv");
	outputFile << "x,y,z,gradX,gradZ" << "\n";

	for (size_t i = 0; i < x.size()-1; i+=2) {
	    	for (size_t j = 0; j < y.size(); j++) {
	         	for (size_t k =0 ; k < z.size()-1; k+=2){
              			outputFile << x[i] << "," << y[j] << "," << z[k] << "," << xGrad[i][j][k] << "," << yGrad[i][j][k] << "\n";
	    	        }
	       	}
	 }
	outputFile.close();
}


//Writes the computed integral to a .csv file
void writeFieldToCSV(const std::vector<double>& x, 
		     const std::vector<double>& y, 
		     const std::vector<double>& z, 
		     std::vector<std::vector<std::vector<double>>>& field){

	std::ofstream outputFile("gradient.csv");
	outputFile << "x,y,z,gradientValue" << "\n";

	for (size_t i = 0; i < x.size()-1; i+=2) {
	    	for (size_t j = 0; j < y.size(); j++) {
	         	for (size_t k =0 ; k < z.size()-1; k+=2){
              			outputFile << x[i] << "," << y[j] << "," << z[k] << "," << field[i][j][k] << "\n";
	    	        }
	       	}
	 }
	outputFile.close();
}

//This will write the position of all deposits to a .csv file so that it may be plotted.

 std::vector<double> arange(double start, double stop, double stepSize){
 	std::vector<double> array;

	 for (double coordinate = start; coordinate <= stop; coordinate += stepSize) {
          	 array.push_back(coordinate);
      	 }
	return array;
 }
