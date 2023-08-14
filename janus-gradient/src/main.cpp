#include <iostream>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <random>
using namespace std;

//Janus particle gradient computation:
std::random_device rd;
std::mt19937 gen(rd());



//Constants:


	const long double particleRadius        = 2*pow(10,-6);

	const long double particleRadiusSquared = pow(particleRadius,2);
	const long double pi	       		= 3.14159265358979323846;
	const long double twoPi			= 2*pi;
	const long double bounds 		= particleRadius*5;
	const long double Intensity		= 200*pow(10,-3); // milliwatt laser
	
	const long double I0			= 2*Intensity/pow(bounds*2,2);
	const long double alphaIron		= 2.43*pow(10,7);  //
	const int 	  nDeposits		= 10;
	const long double kWater 		= 0.598; //W/m·K
	const long double volumePerDeposit      = 4.1666*(10,-4)/nDeposits;
	const long double lambda		= particleRadius;  //wavelength of laser.
	long double stepSize;
	long double dv;
	
	
struct Point{
	long double x;
	long double y;
	long double z;

};


void writeToCSV(vector<long double> x,vector<long double> y,vector<long double> z,vector<vector<vector<long double>>> field);



float fastSqrt(const float n) 
{
   static union{int i; float f;} u;
   u.i = 0x5F375A86 - (*(int*)&n >> 1);
   return (int(3) - n * u.f * u.f) * n * u.f * 0.5f;
}

long double getq(long double x) {
	//Compute 
	long double I = I0 * (1+ cos(  twoPi*x/lambda));
	
	long double output = I/(volumePerDeposit*kWater);
	return output;
	
}	



long double integral(long double x, long double y, long double z,vector<Point> deposits){
	long double contributionSum = 0.0;


    	for (int i = 0; i < deposits.size(); i++){
    	
    		long double distance = sqrt((x-deposits[i].x)*(x-deposits[i].x) + (y-deposits[i].y)*(y-deposits[i].y) + (z-deposits[i].z)*(z-deposits[i].z));
    		
		long double q = getq(deposits[i].x)/distance;
		contributionSum -=  q;
	}

    	return dv*contributionSum;
	
	
}





int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();
	uniform_real_distribution<double> dis(0.0,1.0);
	const long double resolution     = stof(argv[1]);
	stepSize       			 = bounds/(resolution);
	const long double coating        = stof(argv[2])*pi;
	
	dv	 = stepSize*stepSize*stepSize;

	cout<<"step size = " <<stepSize<<endl;


	//Generate coordinates for FeO deposits:
	vector<Point> deposits;
	Point newPoint;
	long double r,theta,phi;
	r = particleRadius;
	for(int i = 0; i < nDeposits; i++){
		theta = dis(gen)*2*pi;
		phi = dis(gen)*pi;
		
		newPoint.x = r*sin(theta)*cos(phi);
		newPoint.y = r*sin(theta)*cos(phi);
		newPoint.z = r*cos(theta);
		deposits.push_back(newPoint);
		
		
	}
	
	
	
	
	
	
	//Initialize space in cartesian coordinates
	 vector<long double> z;
	 for (long double coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
         	z.push_back(coordinate);
     	 }
	 vector<long double> y = {0.0};
	
	 vector<long double> x = z;
	 vector<vector<vector<long double>>> field(x.size(), vector<vector<long double>>(y.size(), vector<long double>(z.size())));
 	 cout<<"Finished initialization"<<endl;
 	
		
		
    		for (int i = 0; i<x.size(); i++){
    			for(int j = 0; j<y.size(); j++){
    				for(int k = 0; k<z.size(); k++){
					//For a given point in 3D space, compute the integral given by Agnese
					field[i][j][k] = integral(x[i],y[j],z[k],deposits);
    				}
    			}
    		}
    	

 
    		
    		

 	
 	
 	

 	
 	float total_points = y.size()*z.size()*x.size();
 	cout<<"total no. of points: "<<total_points<<endl;
 	
	
	//Finite Difference Iterator/////////////
	//Note, this takes alot of time..../////
	/*
	for (int t = 0; t < iterationLimit; t++) {
   	 vector<Point> tempGradientList;



		
    		for ( auto point1 : gradientList) {
    			Point updatedPoint = point1; 
     	 		bool updated = false; 
     	 		
       		  	for ( auto point2 : gradientList) {
            	   		if (point1.x != point2.x && point1.y != point2.y && point1.z != point2.z && point1.r >= particleRadiusSquared) {
            	   		
               				long double dx = (point1.x - point2.x);
               				long double dy = (point1.y - point2.y); 
                			long double dz = (point1.z - point2.z);
               				const long double distance = fastSqrt(dx * dx + dy * dy + dz * dz);
               				
               				
              				if (distance <= (stepSize)) {
                   				updatedPoint.gradientValue +=(point2.gradientValue)/6;
                   				updated = true;
                    				break;
               				}
            			}
       		  	}
     			tempGradientList.push_back(updated ? updatedPoint : point1);
     		}
   		gradientList = move(tempGradientList);
   		cout<<"finished iteration " << t<<endl;
	}*/

	////////////////////////////////////////////////

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
	writeToCSV(x,y,z,field);
		

	///////////////////Compute elapsed time/////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Program completed after: " << elapsed_seconds << " seconds" << std::endl;
	//////////////////////////////////////////////////////////////////S//////////	
	
	return 0;
}




void writeToCSV(vector<long double> x,vector<long double> y,vector<long double> z,vector<vector<vector<long double>>> field) {

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



