#include "functions.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <random>
#include <omp.h>
using namespace std;
std::random_device rd;
std::mt19937 gen(rd());


	

void generateDeposits(vector<Point> &deposits, int nDeposits) {
	size_t depositCounter = 0;
    	uniform_real_distribution<double> phi(0.0,twoPi);
    	uniform_real_distribution<double> costheta(0.0,1);
    	uniform_real_distribution<double> u(0,1);
	
    	while (depositCounter < nDeposits) {
    		long double theta = acos(costheta(gen));
    		long double r 	  = particleRadius* pow(u(gen),-3);
    	
    		long double x = r*sin(theta) * cos(phi(gen));
    		long double y = r*sin(theta) * sin(phi(gen));
    		long double z = r*cos(theta);
   		if ((x*x + y*y + z*z )< particleRadius*particleRadius){
    			deposits.emplace_back(Point{x,y,z});
    		   	depositCounter +=1;
    		}
    		
    		
    	}
}

void writeFieldToCSV(const std::vector<long double>& x, 
		     const std::vector<long double>& y, 
		     const std::vector<long double>& z, 
		     std::vector<std::vector<std::vector<long double>>>& field)
		     {

	std::ofstream outputFile("gradient.csv");
	outputFile << "x,y,z,gradientValue" << "\n";
    	
	for (int i = 0; i < x.size(); i++) {
	    	for (int j = 0; j < y.size(); j++) {
	         	for (int k = 0; k < z.size(); k++){
	                  	outputFile << x[i] << "," << y[j] << "," << z[k] << "," << field[i][j][k] << "\n";
	    	        }
	       	}
	 }
    
	outputFile.close();
}
void writeDepositToCSV(vector<Point> &deposits) {
    std::ofstream outputFile("deposits.csv");
    outputFile << "x,y,z" << "\n";
    
	for (int i = 0; i < size(deposits); i++) {
	        outputFile << deposits[i].x << "," << deposits[i].y  << "," << deposits[i].z  << "\n";
	}

    outputFile.close();
}
