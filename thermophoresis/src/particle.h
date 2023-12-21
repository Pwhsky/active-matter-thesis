/*
Alex Lech 2023

This header contains all declarations for the Particle class, as well as 
other functions used, such as writing to files.
*/
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>
#include <iostream>

struct Point{
	double x;
	double y;
	double z;
};

class Particle{
	public:
		Point  center;
		double radius;
		double selfPropulsion[3];
		std::vector<Point> deposits;
		Particle(Point particleCenter,double r,double velocity){
			center = particleCenter;
			radius = r;
			selfPropulsion[0] = velocity;
			selfPropulsion[1] = velocity;
			selfPropulsion[2] =velocity;
		}  // Constructor with parameters
		//Kinematics
		void updatePosition();
		double getRadialDistance(Point r);
		void rotate(double theta);

		void generateDeposits(int nDeposits);
		void writeDepositToCSV();
		

};


std::vector<double> arange(double start, double stop, double stepSize);
std::vector<Particle> initializeParticles();
void generateConfiguration(std::vector<Point> &deposits,int nDeposits); //Custom configuration (NOT IN USE)

//Writes temperature increase!
void writeFieldToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::vector<std::vector<std::vector<double>>>& field);

//Writes X and Z gradient!
void writeGradToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::vector<std::vector<std::vector<double>>>& xGrad, std::vector<std::vector<std::vector<double>>>& yGrad);



#endif
