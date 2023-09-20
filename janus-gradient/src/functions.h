#include <iostream>
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>

	constexpr long double pi	      		 = 3.14159265358979323846;
	constexpr long double twoPi			 = 2*3.14159265358979323846;
	constexpr long double particleRadius   		 = 2    *pow(10,-6);
	constexpr long double particleRadiusSquared      = particleRadius*particleRadius;
	constexpr long double depositRadius	 	 = 30   *pow(10,-9);	
	constexpr long double volumePerDeposit		 = 4*pi *pow((depositRadius),3)/3; 
	constexpr long double depositArea	 	 = 2*pi *pow(depositRadius,2); 
	

	constexpr long double intensity		  	 = 100  *pow(10,-3);  //milliwatt laser
	constexpr long double areaOfIllumination 	 = 40   *pow(10,-6);  //How much area the laser is distributed on.
	constexpr long double I0		 	 = 2*intensity/(pow(areaOfIllumination*2,2));
	constexpr long double waterConductivity	 = 0.606;
struct Point{
	long double x;
	long double y;
	long double z;
};

void generateDeposits(std::vector<Point> &deposits,int nDeposits);
void writeDepositToCSV(std::vector<Point> &deposits);
void writeFieldToCSV(const std::vector<long double>& x, 
			    const std::vector<long double>& y,
			    const std::vector<long double>& z,
			    std::vector<std::vector<std::vector<long double>>>& field);


#endif
