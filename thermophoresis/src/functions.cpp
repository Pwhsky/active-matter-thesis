#include "functions.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <random>

using namespace std;
	constexpr double particleRadius   		  = 2    *pow(10,-6);


//Places particles next to each other	
std::vector<Particle> initializeParticles(){

	Point centerOfParticle1 = {1.5*particleRadius,0.0,0.0}; 
	Point centerOfParticle2 = {-1.1*particleRadius,0.0,0.0}; 
	Particle particle1(centerOfParticle1,particleRadius,0.0);
	Particle particle2(centerOfParticle2,particleRadius,0.0);
	vector<Particle> particles = {particle1,particle2};

	
	return particles;
	
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

	std::ofstream outputFile("temperature.csv");
	outputFile << "x,y,z,temperature" << "\n";

	for (size_t i = 0; i < x.size()-1; i+=2) {
	    	for (size_t j = 0; j < y.size(); j++) {
	         	for (size_t k =0 ; k < z.size()-1; k+=2){
              			outputFile << x[i] << "," << y[j] << "," << z[k] << "," << field[i][j][k] << "\n";
	    	        }
	       	}
	 }
	outputFile.close();
}

//This will write the position of ALL deposits to a .csv file so that it may be plotted.

 std::vector<double> arange(double start, double stop, double stepSize){
 	std::vector<double> array;

	 for (double coordinate = start; coordinate <= stop; coordinate += stepSize) {
          	 array.push_back(coordinate);
      	 }
	return array;
 }
