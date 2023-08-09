#include <iostream>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>


using namespace std;

//Janus particle gradient computation:


//Constants:

	const long double particleRadius        = pow(10,-3);
	const long double particleRadiusSquared = pow(particleRadius,2);
	const long double pi	       		= 3.14159265358979323846;
	const long double bounds 		= particleRadius*20;
	
	const long double alphaSilica		 = 1*pow(10,-7);  //actual magnitude is 10e-7
	const long double alphaFluid		 = 1.43*pow(10,-7);
	
struct Point {
	long double x, y, z, r;
	long double  gradientValue = 0.0;
};
Point newPoint;
	


long double sumContributions(std::vector<Point> goldCoating, const long double ix, const long double iy, const long double iz);
void writeToCSV(vector<long double> x,vector<long double> y,vector<long double> z,vector<vector<vector<long double>>> field);
void appendValues(Point &newPoint, std::vector<Point>& gradient,long double ix, long double iy, long double iz, long double gradientValue);


float fastSqrt(const float n) 
{
   static union{int i; float f;} u;
   u.i = 0x5F375A86 - (*(int*)&n >> 1);
   return (int(3) - n * u.f * u.f) * n * u.f * 0.5f;
}

long double getq(long double x,long double y, long double z) {
	//The potential will vary depending on the 2 domains:
	//water   r > particleRadius
	//gold    r = particleRadius && z > 0
	//silica  r = particleRadius
	
	
	long double output = exp(-2*particleRadius/pow((r-particleRadius),2)  );

	return output;
	
}	


long double integral(int ix, int iz, int iy, vector<long double> x,vector<long double> y, vector<long double> z,const long double dxdydz){
	long double contributionSum = 0.0;
	long double r 	  =  x[ix]*x[ix]  + z[iz]*z[iz];

	
	//distance = r-r'
	//Volume element = dxdydz
		
    	for (int i = 0; i<x.size(); i++){
    			for(int j = 0; j<y.size(); j++){
    				for(int k = 0; k<z.size(); k++){
    				
				
					long double rPrim = sqrt(x[i]*x[i] + z[k]*z[k]);	
					long double rDiff = r - rPrim;

					long double q = getq(rPrim);
					
					contributionSum -=    ( q*rDiff/(abs(rDiff*rDiff*rDiff)) )*dxdydz;
							
    			}
    		}
    	}
    	return contributionSum/(4*pi);
	
	
}





int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();
	
	const long double resolution     = stof(argv[1]);
	const long double stepSize       = 2*bounds/(resolution);
	const long double coating        = stof(argv[2])*pi;
	const long double dxdydz	 = pow(stepSize,3);

	//Fill list with location of gold coated points and write to csv.

	
	
	vector<Point> gradientList;
	
	
		//Initialize space in spherical coordinates:
	vector<long double> z;
	for (long double coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
        	z.push_back(coordinate);
   
    	}
	vector<long double> y = {0.0};
	vector<long double> x = z;
	



	
	//Initialize space:
	vector<vector<vector<long double>>> field(x.size(), vector<vector<long double>>(y.size(), vector<long double>(z.size())));
 	cout<<"Finished initialization"<<endl;
 	
		//Volume element = r*r*dr*dphi*dtheta 


		
    		for (int i = 0; i<x.size(); i++){
    			for(int j = 0; j<y.size(); j++){
    				for(int k = 0; k<z.size(); k++){
					//For a given point in 3D space, compute the integral given by Agnese
					field[i][j][k] = integral(i,j,k,x,y,z,dxdydz);
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




