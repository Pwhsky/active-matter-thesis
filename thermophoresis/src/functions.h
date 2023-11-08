#include <iostream>
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>

	constexpr double pi	      		 	  	  = 3.14159265358979323846;
	constexpr double twoPi					  = 2*3.14159265358979323846;
	constexpr double particleRadius   		  = 2    *pow(10,-6);
	constexpr double particleRadiusSquared    = particleRadius*particleRadius;
	constexpr double depositRadius	 	      = 30   *pow(10,-9);	
	constexpr double volumePerDeposit	      = 4*pi *pow((depositRadius),3)/3; 
	constexpr double depositArea	 	      = 2*pi *pow(depositRadius,2); 
	constexpr double intensity		  		  = 100  *pow(10,-3);  //Watts
	constexpr double areaOfIllumination 	  = 40   *pow(10,-6);  //Meters  How much area the laser is distributed on.
	constexpr double I0		 				  = 2*intensity/(pow(areaOfIllumination*2,2)); 
	constexpr double waterConductivity	 	  = 0.606;

struct Point{
	double x;
	double y;
	double z;
};

class Particle{
	public:
		Point  center;
		double radius;
		std::vector<Point> deposits;
		Particle(Point particleCenter,double r){
			center = particleCenter;
			radius = r;
		}  // Constructor with parameters
		double getRadialDistance(Point r);
		void generateDeposits(int nDeposits);
		void writeDepositToCSV();
		void rotate(double theta);
};


std::vector<double> arange(double start, double stop, double stepSize);

void generateConfiguration(std::vector<Point> &deposits,int nDeposits); //Custom configuration (NOT IN USE)

//Writes temperature increase!
void writeFieldToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::vector<std::vector<std::vector<double>>>& field);

//Writes X and Z gradient!
void writeGradToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::vector<std::vector<std::vector<double>>>& xGrad, std::vector<std::vector<std::vector<double>>>& yGrad);



#endif
