#include <iostream>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <random>
#include <omp.h>
using namespace std;
std::random_device rd;
std::mt19937 gen(rd());

//Water conductivity = 0.606 

	
	constexpr long double pi	       	 = 3.14159265358979323846;
	constexpr long double twoPi		 = 2*3.14159265358979323846;
	constexpr long double particleRadius     = 2*pow(10,-6);
	constexpr int 	      nDeposits		 = 5000;		
	constexpr long double depositRadius	 = 30*pow(10,-9);	
	
	constexpr long double volumePerDeposit	 = 4*pi*pow((depositRadius),3)/3; 
	constexpr long double depositArea	 = pi*pow(depositRadius,2);
	
	
	constexpr long double bounds 		 = 5*pow(10,-6);     //Area to simulate
	constexpr long double lambda		 = 8000*pow(10,-9);  //wavelength of laser. 
	constexpr long double intensity		 = 100*pow(10,-3);   // milliwatt laser
	constexpr long double areaOfIllumination = 40*pow(10,-6);    //How much area the laser is distributed on.
	constexpr long double I0		 = 2*intensity/(pow(areaOfIllumination*2,2));
	
	long double stepSize;
	long double dv;

	
struct Point{
	long double x;
	long double y;
	long double z;
};

inline void generateDeposits(vector<Point> &deposits);
inline void writeDepositToCSV(vector<Point> &deposits);
inline void writeFieldToCSV(const vector<long double>& x, 
			    const vector<long double>& y,
			    const vector<long double>& z,
			    vector<vector<vector<long double>>>& field);
			    
			    

inline  long double integral(long double x, long double y, long double z,vector<Point> deposits){


	
	long double laserTerm 		        = I0 + I0*cos(twoPi*x/lambda);
	long double qTerm 			= laserTerm*depositArea/(volumePerDeposit);
	
	long double contributionSum = 0.0;
	long double q = 0;

	long double inv_sqrt_distance1;
	//float distance2;
	
    	for (size_t i = 0; i < deposits.size(); i++){
    		inv_sqrt_distance1 = 1.0/sqrt((x-deposits[i].x)*(x-deposits[i].x) + 
    					    (y-deposits[i].y)*(y-deposits[i].y) +
    					    (z-deposits[i].z)*(z-deposits[i].z));
  		
		contributionSum +=  inv_sqrt_distance1;
	}
    	return contributionSum*dv*qTerm/(4*pi*0.606);
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
	 vector<long double> y = {0.0};
	      
	
    const int totalIterations = x.size() * y.size() * z.size();
    int currentIteration = 0;
     	 
	
	 vector<vector<vector<long double>>> field(x.size(), vector<vector<long double>>(y.size(), vector<long double>(z.size())));
 	 cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;


	//As the integral is a triple nested for-loop, parallelization is in order...
 	#pragma omp parallel for
    	for (int i = 0; i<x.size(); i++){
    		for(int j = 0; j<y.size(); j++){
    			for(int k = 0; k<z.size(); k++){

				field[i][j][k] = integral(x[i],y[j],z[k],deposits);
				
               			currentIteration++;

              			  // Calculate progress percentage so that the user has something to look at
              			if(currentIteration % 500 == 0) {
             			  	float progress = round(static_cast<float>(currentIteration) / totalIterations * 100.0);

               				 // Print progress bar
               				 #pragma omp critical
               				 {
                			 		cout << "Progress: "
                			    		<< progress 
                			    		<< "% ("
                			   		<< currentIteration
                			   	 	<< "/" 
                			   	 	<< totalIterations 
                			   		 << ")\r";
                	  	 	 		cout.flush();
                				}
              			}
			}
    		}
    	}
    
	cout<<"Simulation finished, writing to csv..."<<endl;
	writeFieldToCSV(x,y,z,field);
	writeDepositToCSV(deposits);	

	///////////////////Compute elapsed time/////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Program completed after: " << elapsed_seconds << " seconds" << std::endl;
	//////////////////////////////////////////////////////////////////S//////////	
	
	return 0;
}

inline void generateDeposits(vector<Point> &deposits) {
	size_t depositCounter = 0;
    	uniform_real_distribution<double> phi(0.0,twoPi);
    	uniform_real_distribution<double> costheta(0,1);
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
void writeFieldToCSV(const std::vector<long double>& x, const std::vector<long double>& y, const std::vector<long double>& z, std::vector<std::vector<std::vector<long double>>>& field) {

    std::ofstream outputFile("gradient.csv");
    outputFile << "x,y,z,gradientValue" << std::endl;
    
	for (int i = 0; i < x.size(); i++) {
	    	for (int j = 0; j < y.size(); j++) {
	         	for (int k = 0; k < z.size(); k++) {
	                  	outputFile << x[i] << "," << y[j] << "," << z[k] << "," << field[i][j][k] << std::endl;
	    	        	}
	       	}
	 }
    
    outputFile.close();
}

void writeDepositToCSV(vector<Point> &deposits) {

    std::ofstream outputFile("deposits.csv");
    outputFile << "x,y,z" << std::endl;
    
	for (int i = 0; i < nDeposits; i++) {
	        outputFile << deposits[i].x << "," << deposits[i].y  << "," << deposits[i].z  << std::endl;
	}

    outputFile.close();
}


//         const long double volumePerDeposit             = 4.1666*(10,-4)/nDeposits;




/* OLD VERSION, USING FIBONACCI SPHERE.
inline void generateDeposits(vector<Point> &deposits) {
    // Places out the deposits evenly according to the Fibonacci sphere.
    const long double phi     = M_PI * (3.0 - sqrt(5.0)); // Golden ratio constant
    const long double offset = 2.0 / nDeposits; // Offset to distribute points

    for (int i = 0; i < nDeposits; i++) {
    const    long double y = 1 - (i * offset); // y ranges from 1 to -1
    const    long double radius = sqrt(1 - y * y);
        
     const   long double theta = phi * i; // Use golden ratio to distribute points
     const   long double x = cos(theta) * radius * particleRadius;
     const   long double z = sin(theta) * radius * particleRadius;
     
     deposits.emplace_back(Point{x,y*particleRadius,z});
        
    }
}

//New version, sample uniform random numbers untill all are within radius.
*/
