#include "functions.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <random>

using namespace std;
std::random_device rd;
std::mt19937 gen(rd());

	

void Particle::generateDeposits(vector<Point> &deposits, int nDeposits) {
	//Create random number generators, with certain intervals.
	//to create a janus-particle distribution, costheta parameter should be between 0 and 1 (corresponding to z axis).
	//If one wants a completely covered particle, set costheta to (-1,1).
	//To adjust the areas of initialization, play around with the phi and costheta parameters.
	//u is commonly set (0.9,1) so the deposits are near the surface.

    uniform_real_distribution<double> phi(0.0,pi/4); 
    uniform_real_distribution<double> costheta(0.1,1);
    uniform_real_distribution<double> u(0.8,1);

	//Initiate deposits
	for(int i = 0; i<nDeposits; i++){

		double theta = acos(costheta(gen));
		double r 	 = (this->radius)*u(gen);

		//Convert to cartesian:
    	double x = r*sin(theta) * cos(phi(gen)) + this->particleCenter.x; 
    	double y = r*sin(theta) * sin(phi(gen)) + this->particleCenter.y;
    	double z = r*cos(theta)					+ this->particleCenter.z;
   		
		//Add to deposits list
    	deposits.emplace_back(Point{x,y,z});
    	
    }

	//Some interresting configurations:
	//phi(0.0,pi/2)
	//costheta(0.5,0.7)
	//u(0.9,1)

	//phi(0.0,pi/4)
	//costheta(0.5,0.7)
	//u(0.9,1)
}
bool Particle::isOutside(Point r){
	if ((pow(this->particleCenter.x-r.x,2) + 
		pow(this->particleCenter.y-r.y,2) +
		pow(this->particleCenter.z-r.z,2) )> (pow(this->radius,2) )) {
		return true;
	} else{
		return false;
	}

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
void writeDepositToCSV(vector<Point> &deposits) {
    std::ofstream outputFile("deposits.csv");
    outputFile << "x,y,z" << "\n";
    	
	for (size_t i = 0; i < size(deposits); i++) {
	        outputFile << deposits[i].x << "," << deposits[i].y  << "," << deposits[i].z  << "\n";
	}

    outputFile.close();
}
 std::vector<double> arange(double start, double stop, double stepSize){
 	std::vector<double> array;

	 for (double coordinate = start; coordinate <= stop; coordinate += stepSize) {
          	 array.push_back(coordinate);
      	 }
	return array;
 }
