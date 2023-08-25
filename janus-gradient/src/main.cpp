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
	const long double pi	       			       = 3.14159265358979323846;
	const long double twoPi			       = 2*pi;


//Particle:
	const long double particleRadius            	= 2*pow(10,-6);
	const int 	  nDeposits			      		= 1000;			  //should be between 1500-5000 depending on deposit size.
//         const long double volumePerDeposit             = 4.1666*(10,-4)/nDeposits;
	const long double depositRadius			= 10*pow(10,-9);	
	const long double volumePerDeposit		= 4*pi*pow((depositRadius),3)/3; 
	const long double depositArea			= pi*pow(depositRadius,2);
	const long double kSilica				= 1.1;  			//W/m*K
//Simulation box:
	const long double bounds 		   		= 5*pow(10,-6); //40 microns for now
	long double stepSize;
	long double dv;
	const long double alphaWater 			= 1.43*pow(10,-7);

	
//Laser:
	const long double lambda			= 800*pow(10,-9);  //wavelength of laser. 
	const long double intensity			= 100*pow(10,-3); // milliwatt laser
	const long double areaOfIllumination     = 40*pow(10,-6); 
	const long double I0				= 2*intensity/(pow(areaOfIllumination*2,2));
	
	const long double laserTerm 		= I0 + I0*cos(twoPi/lambda);
	const long double qTerm 			= laserTerm*depositArea/(volumePerDeposit);
	
struct Point{
	long double x;
	long double y;
	long double z;
};


inline void writeToCSV(const vector<long double>& x, const vector<long double>& y, const vector<long double>& z, vector<vector<vector<long double>>>& field);

float fastSqrt(const float n);


inline long double integral(long double x, long double y, long double z,vector<Point> deposits){
	float contributionSum = 0.0;
	float q = 0;
	//long double positionFromCenter = sqrt(x*x + y*y + z*z);
	float distance1;
	//float distance2;
	int i;

    	for (i = 0; i < deposits.size(); i++){
    		distance1 = sqrt((x-deposits[i].x)*(x-deposits[i].x) + 
    					    (y-deposits[i].y)*(y-deposits[i].y) +
    					    (z-deposits[i].z)*(z-deposits[i].z));
  		

		
		contributionSum +=  1/distance1;

	}
    	return contributionSum*dv*qTerm;
}


inline void generateDeposits(vector<Point> &deposits) {
    // Places out the deposits evenly according to the Fibonacci sphere.
    Point newPoint;
    const long double phi = M_PI * (3.0 - sqrt(5.0)); // Golden ratio constant
    const long double offset = 2.0 / nDeposits; // Offset to distribute points

    for (int i = 0; i < nDeposits; i++) {
    const    long double y = 1 - (i * offset); // y ranges from 1 to -1
    const    long double radius = sqrt(1 - y * y);
        
     const   long double theta = phi * i; // Use golden ratio to distribute points
     const   long double x = cos(theta) * radius * particleRadius;
     const   long double z = sin(theta) * radius * particleRadius;

        newPoint.x = x;
        newPoint.y = y * particleRadius;
        newPoint.z = z;
        deposits.push_back(newPoint);
        
    }
}



int main(int argc, char** argv) {
	auto startTimer = std::chrono::high_resolution_clock::now();

	const long double resolution     = stof(argv[1]);
	stepSize       			     = bounds/(resolution);
					 
	
	dv = stepSize*stepSize*stepSize;
	cout<<"step size = " <<stepSize<<endl;

//Generate coordinates for FeO deposits:
	vector<Point> deposits;
	generateDeposits(deposits);
	
	
//Initialize space in cartesian coordinates
	vector< long double> spaceVector;
	for (long double coordinate = -bounds; coordinate <= bounds; coordinate += stepSize) {
        	spaceVector.push_back(coordinate);
     	}
     	const vector<long double> z = spaceVector;
	const vector<long double> x = spaceVector;
	const vector<long double> y = {0.0};
	      
	if (stof(argv[2]) == 3) {
		const vector<long double> y = z;
	}

     	 
    const int totalIterations = x.size() * y.size() * z.size();
    int currentIteration = 0;
     	 
	
	 vector<vector<vector<long double>>> field(x.size(), vector<vector<long double>>(y.size(), vector<long double>(z.size())));
 	 cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
 	
 	#pragma omp parallel for
    	for (int i = 0; i<x.size(); i++){
    		for(int j = 0; j<y.size(); j++){
    			for(int k = 0; k<z.size(); k++){
    			
    			
				//For a given point in 3D space, compute the integral given by Agnese
				field[i][j][k] = integral(x[i],y[j],z[k],deposits);
				
				                // Increment the current Iteration counter
                currentIteration++;

                // Calculate progress percentage
                if(currentIteration % 500 == 0) {
               			
             		 float progress = round(static_cast<float>(currentIteration) / totalIterations * 100.0);

               	 // Print progress bar
               		 #pragma omp critical
               	 			{
                			    std::cout << "Progress: " << progress << "% (" << currentIteration << "/" << totalIterations << ")\r";
                	  	 	 std::cout.flush();
                			}
              		  	}
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

/*
float fastSqrt(const float n) 
{
   static union{int i; float f;} u;
   u.i = 0x5F375A86 - (*(int*)&n >> 1);
   return (int(3) - n * u.f * u.f) * n * u.f * 0.5f;
}

*/
void writeToCSV(const std::vector<long double>& x, const std::vector<long double>& y, const std::vector<long double>& z, std::vector<std::vector<std::vector<long double>>>& field) {
    std::ofstream outputFile("gradient.csv");
    outputFile << "x,y,z,gradientValue" << std::endl;
    /*
    std::string buffer;

    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < y.size(); j++) {
            for (int k = 0; k < z.size(); k++) {
                buffer += std::to_string(x[i]) + "," + std::to_string(y[j]) + ","
                          + std::to_string(z[k]) + "," + std::to_string(field[i][j][k]) + "\n";
            }
        }
    }

    // Write the entire buffer to the file in a single operation
    outputFile << buffer;
*/
   
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < y.size(); j++) {
            for (int k = 0; k < z.size(); k++) {
                    outputFile << x[i] << "," << y[j] << "," << z[k] << "," << field[i][j][k] << std::endl;
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


