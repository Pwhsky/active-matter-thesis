#include <iostream>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <random>
#include <omp.h>
using namespace std;

//Janus particle gradient computation:
std::random_device rd;
std::mt19937 gen(rd());



//Constants:
	const long double pi	       		= 3.14159265358979323846;
	const long double twoPi			= 2*pi;


	//Particle:
	const long double particleRadius        = 2*pow(10,-6);
	const long double particleRadiusSquared = pow(particleRadius,2);
	const int 	  nDeposits		= 6500;
	const long double volumePerDeposit      = 4.1666*(10,-4)/nDeposits;
	//const long double volumePerDeposit	= 4*pi*pow((100*pow(10,-9)),3)/3;
	
	
	//Simulation box:
	const long double bounds 		= particleRadius*20; //40 microns for now
	const long double kWater 		= 0.598; //W/m·K
	      long double dv;
	      long double stepSize;
	      
	//Laser:
	const long double lambda		= 808*pow(10,-9);  //wavelength of laser. Keep above particleRadius generally
	const long double Intensity		= 100*pow(10,-3); // milliwatt laser
	const long double I0			= 2*Intensity/pow(bounds*2,2);

	
	
struct Point{
	long double x;
	long double y;
	long double z;

};


void writeToCSV(vector<long double> x,vector<long double> y,vector<long double> z,vector<vector<vector<long double>>> field);
float fastSqrt(const float n);


inline long double integral(long double x, long double y, long double z,vector<Point> deposits){
	long double contributionSum = 0.0;
	long double q;
    	for (int i = 0; i < deposits.size(); i++){
    		long double distance = sqrt((x-deposits[i].x)*(x-deposits[i].x) + 
    					    (y-deposits[i].y)*(y-deposits[i].y) +
    					    (z-deposits[i].z)*(z-deposits[i].z));
    		
    		
    		
    		
    		q = ((I0 + I0*cos(twoPi*distance/lambda))/volumePerDeposit)/(distance*kWater);
		contributionSum -=  q;
	}
    	return contributionSum*dv;
}


inline void generateDeposits(vector<Point> &deposits) {

//Places out the deposits evenly according to fibonacci sphere.
	Point newPoint;
	long double r,theta,phi;
	phi = pi*(3.0 -sqrt(5.));
	r = particleRadius;
	for(int i = 0; i < nDeposits; i++){
		newPoint.y = (1- pow((i/(nDeposits-1)),2)) *r;
		long double tempRad = sqrt(1-newPoint.y*newPoint.y);
		theta = phi*i;
		newPoint.x = cos(theta) * tempRad*r;
		newPoint.z = sin(theta) * tempRad*r;
		
		deposits.push_back(newPoint);	
	}
}



int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();

	const long double resolution     = stof(argv[1]);
	stepSize       			 = bounds/(resolution);
	
	
	dv = stepSize*stepSize*stepSize;
	cout<<"step size = " <<stepSize<<endl;

//Generate coordinates for FeO deposits:
	vector<Point> deposits;
	generateDeposits(deposits);
	
	
//Initialize space in cartesian coordinates
	vector<long double> spaceVector;
	for (long double coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
        	spaceVector.push_back(coordinate);
     	}
     	const vector<long double> z = spaceVector;
	const vector<long double> x = spaceVector;
	const vector<long double> y = {0.0};
     	 
     	 
	
	 vector<vector<vector<long double>>> field(x.size(), vector<vector<long double>>(y.size(), vector<long double>(z.size())));
 	 cout<<"Finished initialization"<<endl;
 	
 	#pragma omp parallel for
    	for (int i = 0; i<x.size(); i++){
    		for(int j = 0; j<y.size(); j++){
    			for(int k = 0; k<z.size(); k++){
				//For a given point in 3D space, compute the integral given by Agnese
				field[i][j][k] = integral(x[i],y[j],z[k],deposits);
    			}
    		}
    	}
    
	
	cout<<"Simulation finished, writing to csv..."<<endl;
	writeToCSV(x,y,z,field);
		

	///////////////////Compute elapsed time/////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Program completed after: " << elapsed_seconds << " seconds" << std::endl;
	//////////////////////////////////////////////////////////////////S//////////	
	
	return 0;
}


float fastSqrt(const float n) 
{
   static union{int i; float f;} u;
   u.i = 0x5F375A86 - (*(int*)&n >> 1);
   return (int(3) - n * u.f * u.f) * n * u.f * 0.5f;
}


void writeToCSV(const vector<long double> x,const vector<long double> y, const vector<long double> z,vector<vector<vector<long double>>> field) {

	ofstream outputFile("gradient.csv");
   	outputFile << "x,y,z,gradientValue" << std::endl;

   for (int i = 0; i <x.size(); i++){
	for (int j = 0; j < y.size(); j++){
			for (int k = 0; k < z.size(); k++){
				outputFile << x[i] << "," <<y[j] <<"," << z[k]<<"," << field[i][j][k] << std::endl;		
			}
		}
    }
    	outputFile.close();

}

/*
inline long double getq(long double x) {
	//Compute 
	//long double I = I0 * (1+ cos(  twoPi*x/lambda));
	
	
	long double output = (I0 + I0*cos(twoPi*abs(x)/lambda))/volumePerDeposit/kWater;
	return output;
}*/

/*
inline void generateDeposits(vector<Point> &deposits) {

	
	Point newPoint;
	long double r,theta,phi;
	
	for(int i = 0; i < nDeposits; i++){
		theta = dis(gen)*2*pi;
		phi   = dis(gen)*pi;
		r = particleRadius*dis(gen);
		newPoint.x = r*sin(theta)*cos(phi);
		newPoint.y = r*sin(theta)*sin(phi);
		newPoint.z = r*cos(theta);
		deposits.push_back(newPoint);	
	}
}*/


