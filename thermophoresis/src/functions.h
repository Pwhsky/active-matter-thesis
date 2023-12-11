#include <iostream>
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>


struct Point{
	double x;
	double y;
	double z;
};

class Particle{
	public:
		Point  center;
		double radius;
		double selfPropulsion[2];
		std::vector<Point> deposits;
		Particle(Point particleCenter,double r,double velocity){
			center = particleCenter;
			radius = r;
			selfPropulsion[0] = velocity;
			selfPropulsion[1] = velocity;
		}  // Constructor with parameters
		double getRadialDistance(Point r);
		void generateDeposits(int nDeposits);
		void writeDepositToCSV();
		void rotate(double theta);

};


std::vector<double> arange(double start, double stop, double stepSize);
std::vector<Particle> initializeParticles();
void generateConfiguration(std::vector<Point> &deposits,int nDeposits); //Custom configuration (NOT IN USE)

//Writes temperature increase!
void writeFieldToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::vector<std::vector<std::vector<double>>>& field);

//Writes X and Z gradient!
void writeGradToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::vector<std::vector<std::vector<double>>>& xGrad, std::vector<std::vector<std::vector<double>>>& yGrad);



#endif
