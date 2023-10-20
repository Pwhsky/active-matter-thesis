#include "functions.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>
#include <iomanip>
#include <random>
#include <omp.h>
using namespace std;
std::random_device rd;
std::mt19937 gen(rd());

std::string to_string_with_precision(double value, int precision = 15) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}
	

void generateDeposits(vector<Point> &deposits, int nDeposits) {
	size_t depositCounter = 0;
    	uniform_real_distribution<double> phi(0.0,twoPi);
    	uniform_real_distribution<double> costheta(-0.05,1);
    	uniform_real_distribution<double> u(0.95,1);
	
    	while (depositCounter < nDeposits) {
    		double theta = acos(costheta(gen));
    		double r 	  = particleRadius* u(gen);
    	
    		double x = r*sin(theta) * cos(phi(gen));
    		double y = r*sin(theta) * sin(phi(gen));
    		double z = r*cos(theta);
   		if ((x*x + y*y + z*z )< particleRadius*particleRadius){
    			deposits.emplace_back(Point{x,y,z});
    		   	depositCounter +=1;
    		}
    		
    		
    	}
}
void writeFieldToCSV(const std::vector<double>& x, 
		     const std::vector<double>& y, 
		     const std::vector<double>& z, 
		     std::vector<std::vector<std::vector<double>>>& field)
		     {

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

void writeDepositToCSV(vector<Point> &deposits) {
    std::ofstream outputFile("deposits.csv");
    outputFile << "x,y,z" << "\n";
    	
	for (size_t i = 0; i < size(deposits); i++) {
	        outputFile << deposits[i].x << "," << deposits[i].y  << "," << deposits[i].z  << "\n";
	}

    outputFile.close();
}
