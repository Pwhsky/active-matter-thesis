#include <iostream>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

#include <fstream>


using namespace std;

//Janus particle gradient computation:


//Constants:
	const long double scale 	    = pow(10,-3);
	const long double particleRadius    = 1e-2;
	const long double particleRadiusSquared = pow(particleRadius,2);
	const long double pi	            = 3.14159265358979323846;
	//const long double ix		    = 0.0;
	const double viscosity              = 1*pow(10,-3);
	const double stokesCoefficient	    = 6.0*pi*viscosity*particleRadius;			   //Phoretic strength
	const long double bounds 	    = 0.03;

struct Point {
	long double x;
	long double y;
	long double z;
	long double gradientValue;
};


long double sumContributions(std::vector<Point> goldCoating, long double ix, long double iy, long double iz);
void writeToCSV(std::vector<Point> points);
void appendValues(std::vector<Point>& gradient,long double ix, long double iy, long double iz, long double gradientValue);

std::vector<Point> generateCap(long double capResolution,long double coating){

	vector<Point> goldCoating;
	ofstream capFile("cap.csv");
   	capFile << "x,y,z" << std::endl;

   	Point newPoint;
	for (int i = 0; i < capResolution; i++){
		long double theta = i*2*pi/capResolution;
		for (int j = 0; j < capResolution; j++){
			long double phi = j*pi/capResolution;
			
			newPoint.x = particleRadius*cos(theta)*sin(phi);
			newPoint.y = particleRadius*sin(theta)*sin(phi);
			
			if (phi < coating ){ //This controls how much of the sphere is coated. pihalf = half coated.
				newPoint.z = particleRadius*cos(phi);
				goldCoating.push_back(newPoint);
				capFile << newPoint.x << "," << newPoint.y << "," << newPoint.z << std::endl;
			}

   	
		}
	}
	
	
	return goldCoating;
}


int main(int argc, char** argv) {

	const long double resolution  = stof(argv[1]);
	const long double resolutionSize = 1/resolution;
	auto startTimer = std::chrono::high_resolution_clock::now();
	const long double capResolution = 50;
	const long double coating 	  = stof(argv[2])*pi;
	
	//Fill list with location of gold coated points and write to csv.
	vector<Point> goldCoating = generateCap(capResolution,coating);


	
	vector<Point> gradientList;
	for (long double ix = -bounds; ix<bounds; ix += resolutionSize){
		for (long double iy = -bounds; iy<bounds; iy += resolutionSize){
			for (long double iz = -bounds; iz<bounds; iz+= resolutionSize){
		
				long double distance = sqrt(ix*ix+iy*iy+iz*iz); //Distance from origo
					
			//Check if we are inside the particle
				if(distance> particleRadius){
			
					long double contributionSum = sumContributions(goldCoating,ix,iy,iz);
					long double gradientValue = 1/(4*pi*contributionSum);		
					appendValues(gradientList,ix,iy,iz,gradientValue);
				
				//continue;
				
				}
			
				else {appendValues(gradientList,ix,iy,iz,0.0);}
			
			

			}
		}
 	}

		cout<<"Simulation finished, writing to csv..."<<endl;
		writeToCSV(gradientList);
		

	///////////////////Compute elapsed time/////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Program completed after: " << elapsed_seconds << " seconds" << std::endl;
	////////////////////////////////////////////////////////////////////////////	
	
	return 0;
}

long double sumContributions(std::vector<Point> goldCoating, long double ix, long double iy, long double iz) {
	long double contributionSum = 0.0;
	for (auto point : goldCoating){
		long double dx = (point.x-ix);
		long double dy = (point.y-iy);
		long double dz = (point.z-iz);
		long double goldDistance = sqrt(dx*dx + dy*dy  +dz*dz); 
		contributionSum += goldDistance/goldCoating.size();	
	}
	return contributionSum;
}

void appendValues(std::vector<Point>& gradient,long double ix, long double iy, long double iz, long double gradientValue) {
	Point newPoint;
	newPoint.x = ix;
	newPoint.y = iy;
	newPoint.z = iz;
	newPoint.gradientValue = gradientValue;
	gradient.push_back(newPoint);
}

void writeToCSV(std::vector<Point> points) {

	ofstream outputFile("gradient.csv");
	
   	outputFile << "x,y,z,gradientValue" << std::endl;

   	
   	
    // Write field values
	for (auto point : points) {
  		outputFile << point.x << "," << point.y << "," << point.z << "," << point.gradientValue << std::endl;
   	}
    	outputFile.close();
}

