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
		double velocity[3];
		double selfRotation[3];

		std::vector<Point> deposits;
		Particle(Point particleCenter,double r,double vel){
			center = particleCenter;
			radius = r;
			for(int i = 0; i<3;i++){
				velocity[i] = vel;
				selfRotation[i]   = vel;
			}

		}  
		//Kinematics
		double getRadialDistance(Point r);
		void   updatePosition();
		void   rotate(double theta);
		void   rotation_transform();
		void   generateDeposits(int nDeposits);
		void   writeDepositToCSV();
		void   hard_sphere_correction();
};

/*
Example: to create a particle with velocity 0.0 at origin with radius 1 micron, you do:

	Point center = {0.0,0.0,0.0};

	Particle myParticle(center, 1e-6 , 0.0);
*/



std::vector<Particle> initializeParticles();std::vector<double> get_new_coordinates(std::vector<double> omega, std::vector<double> x);



//Linear algebra
std::vector<std::vector<double>> mat_mat_mul (std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);
std::vector<double> mat_vec_mul(std::vector<std::vector<double>> R, std::vector<double> x);
std::vector<double> arange(double start, double stop, double stepSize);
std::vector<double> cross_product(std::vector<double> a,std::vector<double> b);

double get_norm(std::vector<double> a);


void writeFieldToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, 
							std::vector<std::vector<std::vector<double>>>& field);

void writeGradToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, 
							std::vector<std::vector<std::vector<double>>>& xGrad, std::vector<std::vector<std::vector<double>>>& yGrad);

void generateConfiguration(std::vector<Point> &deposits,int nDeposits); //Custom configuration (NOT IN USE)
#endif
