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
	const long double bounds 	    = 0.025;
	const long double capResolution     = 50;


struct Point {
	long double x, y, z, gradientValue;
};
	Point newPoint;
	


long double sumContributions(std::vector<Point> goldCoating, const long double ix, const long double iy, const long double iz);
void writeToCSV(std::vector<Point> points);
void appendValues(Point &newPoint, std::vector<Point>& gradient,long double ix, long double iy, long double iz, long double gradientValue);


float fastSqrt(const float n) 
{
   static union{int i; float f;} u;
   u.i = 0x5F375A86 - (*(int*)&n >> 1);
   return (int(3) - n * u.f * u.f) * n * u.f * 0.5f;
}

/*
std::vector<Point> generateCap(Point &newPoint, const long double capResolution, const long double coating){

	vector<Point> goldCoating;
	ofstream capFile("cap.csv");
   	capFile << "x,y,z" << std::endl;

   	
	for (int i = 0; i < capResolution; i++){
		const long double theta = i*2*pi/capResolution;
		for (int j = 0; j < capResolution; j++){
		
			const long double phi = j*pi/capResolution;
			
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
}*/

int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();
	
	const long double resolution  = stof(argv[1]);
	const long double stepSize    = 1/resolution;
	const long double coating     = stof(argv[2])*pi;
	
	//Fill list with location of gold coated points and write to csv.

	
	
	vector<Point> gradientList;
	
	//Initialize space:
	vector<long double> y;
	for (long double coordinate = 0.0; coordinate <= bounds; coordinate += stepSize) {
        	y.push_back(coordinate);
    	}
	
	vector<long double> x;
    	for (long double coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
        	x.push_back(coordinate);
    	}
	vector<long double> z = x;
	//Initialize gradient values:
	
	
	//Initialize boundary conditions (spherical cap)
	for (const auto iy : y){
		for (const auto ix : x){
			for (const auto iz : z){
				const long double distanceSquared = ix*ix+iy*iy+iz*iz; //Distance from origo	
				
				if(distanceSquared >= particleRadiusSquared-stepSize/8 && distanceSquared <= particleRadiusSquared+stepSize/8)   {
					newPoint.x = ix;
					newPoint.y = iy;
					if(iz>0){
						newPoint.z = iz;
						newPoint.gradientValue = 0.2;
						gradientList.push_back(newPoint);
						}
					
				}else if (distanceSquared >= particleRadiusSquared){
					newPoint.x = ix;
					newPoint.y = iy;
					newPoint.z = iz;
					newPoint.gradientValue = 0.0;
					gradientList.push_back(newPoint);
				
				}
			}
		}
 	}
 	
 	
 	
 	
 	
 	
 	
 	
 	float total_points = y.size()*z.size()*x.size();
 	cout<<"total no. of points: "<<total_points<<endl;
 	
	
	//iterate over time:
	const int iterationLimit = 2;
	for(int t = 0; t<iterationLimit; t++){
		
		//Find adjacent points
 	for (auto& point1 : gradientList){
 		for (auto& point2 : gradientList){
 		
 			if(&point1 != &point2) {
 				long double dx = (point1.x - point2.x);
 				long double dy = (point1.y - point2.y); 
 				long double dz = (point1.z - point2.z);
				const long double distanceSquared = dx*dx+dy*dy+dz*dz;
				
				if (distanceSquared<=(stepSize*stepSize)){
					
					point1.gradientValue = (point2.gradientValue)/6;
				}
				
 			}
 		
 		
 		}
 	}
		
 	}
		

	
	
	
	



		//The optimizer loves const datatypes, sorry for the clutter...
		
		/*
		
	for (long double y = 0.0; y<bounds; y += stepSize){
		const long double iy = y;
		for (long double x = -bounds; x<bounds; x += stepSize){
			const long double ix = x;
			for (long double z = -bounds; z<bounds; z += stepSize){
				const long double iz = z;
				
				const long double distanceSquared = ix*ix+iy*iy+iz*iz; //Distance from origo

				//Check if the point is inside the particle:
				
				if(distanceSquared>particleRadiusSquared){
				
					const long double sum = sumContributions(goldCoating,ix,iy,iz);
					const long double gradientValue = 1/(sum*sum);		
					appendValues(newPoint,gradientList,ix,iy,iz,gradientValue);
				}
			
				
			}
		}
 	}
	*/
	
	cout<<"Simulation finished, writing to csv..."<<endl;
	writeToCSV(gradientList);
		

	///////////////////Compute elapsed time/////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Program completed after: " << elapsed_seconds << " seconds" << std::endl;
	//////////////////////////////////////////////////////////////////S//////////	
	
	return 0;
}

long double sumContributions(std::vector<Point> goldCoating, const long double ix, const long double iy, const long double iz) {

	long double contributionSum = 0.0;
	for (auto point : goldCoating){
		const long double dx = (point.x-ix);
		const long double dy = (point.y-iy);
		const long double dz = (point.z-iz);
		const long double goldDistance = fastSqrt(dx*dx + dy*dy  +dz*dz); 
		contributionSum += goldDistance;	
	}
	return (contributionSum);
}

void appendValues(Point &newPoint,std::vector<Point>& gradient,const long double ix, const long double iy, const long double iz, const long double gradientValue) {
	
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





