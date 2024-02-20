#include "particle.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <random>

using namespace std;
	constexpr double particleRadius	= 2*pow(10,-6);


//Places particles next to each other	
std::vector<Particle> initializeParticles(){

	//Create the coordinates
	Point centerOfParticle1 = {1.75*particleRadius,0.0,0.0}; 
	Point centerOfParticle2 = {-1.75*particleRadius,0.0,0.0}; 
	Point centerOfParticle3 = {0.0,0.0,-1.75*particleRadius}; 
	Point centerOfParticle4 = {0.0,0.0,1.75*particleRadius}; 

	//Create the particles
	Particle particle1(centerOfParticle1,particleRadius,0.0);
	Particle particle2(centerOfParticle2,particleRadius,0.0);
	Particle particle3(centerOfParticle3,particleRadius,0.0);
	Particle particle4(centerOfParticle4,particleRadius,0.0);
	vector<Particle> particles = {particle1,particle2};

	return particles;
	
}

//Writes the computed integral to a .csv file
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


double get_norm(std::vector<double> a){
	return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}


 std::vector<double> arange(double start, double stop, double stepSize){
 	std::vector<double> array;

	 for (double coordinate = start; coordinate <= stop; coordinate += stepSize) {
          	 array.push_back(coordinate);
      	 }
	return array;
 }
 std::vector<double> cross_product(std::vector<double> a,std::vector<double> b){
	std::vector<double> output = {a[1]*b[2]-a[2]*b[1],
								  a[2]*b[0]-a[0]*b[2],
								  a[0]*b[1]-a[1]*b[0]};
	
	return output;
 }
 std::vector<std::vector<double>> mat_mat_mul (std::vector<std::vector<double>> a,
 										  std::vector<std::vector<double>> b){

	std::vector<std::vector<double>> output = a;

	for(int i = 0; i<3;i++){
		for(int j = 0; j<3; j++) {
			double sum = 0.0;
			for(int k = 0; k<3; k++){
				sum += a[i][k]*b[k][j];
			}
				output[i][j] = sum;
		}
	}
	return output;					
}
std::vector<double> mat_vec_mul(std::vector<std::vector<double>> R, std::vector<double> x){
    std::vector<double> temp(3);
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += R[i][j] * x[j]; // Corrected the indexing here
        }
        temp[i] = sum;
    }
    return temp;
}