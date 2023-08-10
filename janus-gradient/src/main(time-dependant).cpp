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
	const long double bounds 	    = 0.015;
	
	const long double alphaSilica	 = 1*pow(10,-7);  //actual magnitude is 10e-7
	const long double alphaFluid	 = 1.43*pow(10,-7);
struct Point {
	long double x, y, z,r;
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



int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();
	
	const long double resolution     = stof(argv[1]);
	const long double stepSize       = 1/(resolution*1.1);
	const long double coating        = stof(argv[2])*pi;
	const int         iterationLimit = stoi(argv[3]);
	

	const long double skinThickness  = stepSize/(16);
	const long double delta 	 = 1/(stepSize*stepSize) * 0.01;
	//Fill list with location of gold coated points and write to csv.

	
	
	vector<Point> gradientList;
	
	
	
	
	//Initialize space:
	
	vector<long double> y;
	for (long double coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
        	y.push_back(coordinate);
    	}
	vector<long double> z = y;
	vector<long double> x = y;
	int xSize = x.size();
	int ySize = y.size();
	int zSize = z.size();
	
	
	vector<vector<vector<long double>>> field(xSize, vector<vector<long double>>(ySize, vector<long double>(ySize)));
 	long double gradientValue = 0.0;
 	//for (int i = 0; i < x.size(); i++){

 	for (int i = 0; i<xSize; i++){
		for (int j = 0; j < ySize; j++){
			for (int k = 0; k < zSize; k++){
				long double ix = x[i];
				long double jy = y[j];
				long double kz = z[k];
				const long double distanceSquared = ix*ix+jy*jy+kz*kz; //Distance from origo	
				
				if (distanceSquared <= particleRadiusSquared && kz>0){
					gradientValue = 0.5;

				}else if (distanceSquared <= particleRadiusSquared) {
					gradientValue = 0.3;
					
				}else {
					gradientValue = 0.0;
				}	
				field[i][j][k] = gradientValue;		
			}
		}
 	}
 	cout<<"Finished initialization"<<endl;
 	

	for(int t = 0; t<iterationLimit; t++){	
		vector<vector<vector<long double>>> tempField(xSize, vector<vector<long double>>(ySize, vector<long double>(ySize)));
    		for (int i = 1; i<(xSize-1); i++){
    			for(int j = 1; j<(ySize-1); j++){
    				for(int k = 1; k<(zSize-1); k++){

    					const long double jy = y[j];
					const long double ix = x[i];
					const long double kz = z[k];
					const long double distanceSquared = ix*ix+jy*jy+kz*kz;
				
				// 
				//Try with second derivative: 	(Works)
					//if(i-1 == 0 || i+1 == xSize-1 || j-1 == 0 || k-1 == 0 || j+1 == ySize-1 ||k+1 == zSize-1){	
					//	tempField[i][j][k] = 0.0; //This affects the cooling

					//}else{
						tempField[i][j][k] *= (i != 0 && i != xSize-1 && j != 0 && k != 0 && j != ySize-1 && k != zSize-1);
						long double gradX = (field[i+1][j][k] - 2*field[i][j][k] + field[i-1][j][k]);
						long double gradY = (field[i][j+1][k] - 2*field[i][j][k] + field[i][j-1][k]);
						long double gradZ = (field[i][j][k+1] - 2*field[i][j][k] + field[i][j][k-1]);
						
						
						
						if (distanceSquared <= particleRadiusSquared) {
							//Use silica heat coefficient when inside the sphere
							tempField[i][j][k]  =  field[i][j][k] + alphaSilica*(gradY+gradZ+gradX)*delta;
						}else{
							//use water lutidine heat constants.
							tempField[i][j][k]  =  field[i][j][k] + alphaFluid*(gradY+gradZ+gradX)*delta;
						}
					
    					//}
    				}
    			}
    		}
    		field = move(tempField);

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




