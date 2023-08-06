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
	const long double coating 	    = 0.5*pi;

	const double kb      		    = 1.380649 * pow(10,-23);
	const double viscosity              = 1*pow(10,-3);
	const double stokesCoefficient	    = 6.0*pi*viscosity*particleRadius;			   //Phoretic strength
	const long double conductivity	    = 0.606; //W/(m*K)
	const long double bounds 	    = 0.03;


struct Point {
	long double x;
	long double y;
	long double z;
	long double gradientValue;
};



void writeToCSV(std::vector<Point> points);


int main(int argc, char** argv) {
	long double resolution  = stof(argv[1]);
	long double resolutionSize = 1/resolution;
	auto startTimer = std::chrono::high_resolution_clock::now();
	long double capResolution = 100;

	
	//Fill list with location of gold coated points and write to csv.
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
	capFile.close();
	
		
		//Sweep across space:
		vector<Point> gradient;
		
		//for (long double ix = -bounds; ix<bounds; ix += resolutionSize){
			long double ix=0.0;
			for (long double iy = -bounds; iy<bounds; iy += resolutionSize){
				for (long double iz = -bounds; iz<bounds; iz+= resolutionSize){
					
					long double distance = ix*ix + iy*iy+iz*iz ;
					
				//Check if we are inside the particle
					if(distance< particleRadiusSquared){
					
						
						newPoint.x = ix;
						newPoint.y = iy;
						newPoint.z = iz;
						newPoint.gradientValue = 0;
						gradient.push_back(newPoint);
	
					}
					else{
					
						long double contributionSum = 0;
						//Sum contributions:
						for (auto point : goldCoating){
							long double dx = (point.x-ix);
							long double dy = (point.y-iy);
							long double dz = (point.z-iz);
							long double goldDistance = sqrt(dx*dx + dy*dy  +dz*dz); 
							contributionSum += goldDistance/goldCoating.size();	
						}
					
						newPoint.x = ix;
						newPoint.y = iy;
						newPoint.z = iz;
						newPoint.gradientValue = 1/(contributionSum*contributionSum);
						gradient.push_back(newPoint);

					} 
				
			
				}
			}	
		//}
		cout<<"Simulation finished, writing to csv..."<<endl;
		writeToCSV(gradient);
		
	
		

		
		
	///////////////////Compute elapsed time/////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Program completed after: " << elapsed_seconds << " seconds" << std::endl;
	////////////////////////////////////////////////////////////////////////////	
	
	return 0;
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

