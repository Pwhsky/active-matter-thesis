/*
Alex Lech 2023	

This source code contains the time evolution for simulating thermophoresis.
*/

#include <chrono>
#include <fstream>
#include <cmath>
#include <vector>
#include "particle.h"


using namespace std;
	int     number_of_steps,nDeposits;
 	double  bounds, lambda, stepSize, dv, dl; 
	bool    onlyTangential = false;
	vector<Point> globalDeposits;


int main(int argc, char** argv) {
	   
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    		= stold(argv[2])  * pow(10,-6); 
	stepSize  		= bounds/(double)(300);		  //Step size, based off of bounds parameter Set to 100 for now
	nDeposits 		= stof(argv[1]);				  //number of deposits to initialize
	lambda	 	    = stold(argv[3])  * pow(10,-9); //Spatial periodicity
	number_of_steps = (int)stof(argv[4]);
    dv	      		= stepSize*stepSize*stepSize;  //volume element for integral

	
	std::ofstream p1("particle_1.csv");
	std::ofstream p2("particle_2.csv");


	p1   << "x,y,z"<<"\n";
	p2   << "x,y,z"<<"\n";

	vector<Particle> particles = initializeParticles();

	for(int i = 0; i < particles.size(); i++ ){
		particles[i].generateDeposits(nDeposits);
		for(int j = 0; j<nDeposits;j++){
			globalDeposits.push_back(particles[i].deposits[j]);
		}
	} 


	//Generate linspace vectors:
	vector<double> linspace         = arange(-bounds,bounds,stepSize);
 	cout<<"Finished initialization of "<< nDeposits <<" deposits."<<endl;
	dl = stepSize*stepSize;
	double thickness = pow(10*stepSize,2);

	
	//Write initial positions:
	for(int i = 0; i<particles.size(); i++){

	}
	p1 << particles[0].center.x << "," << particles[0].center.y << "," << particles[0].center.z << "\n";
	p2 << particles[1].center.x << "," << particles[1].center.y << "," << particles[1].center.z << "\n";

	

	for(int time = 0; time < number_of_steps; time ++){ 

		for(auto &particle:particles){
			particle.getKinematics(linspace,thickness,dl,globalDeposits,lambda,dv);
			particle.rotation_transform();
		}

		for(auto &particle:particles){
			particle.update_position();
			particle.brownian_noise();
		}
		hard_sphere_correction(particles);
		globalDeposits = update_globalDeposits(particles);


			p1 << particles[0].center.x << "," << particles[0].center.y << "," << particles[0].center.z << "\n";
			p2 << particles[1].center.x << "," << particles[1].center.y << "," << particles[1].center.z << "\n";
			//TODO: make exporting particle positions a trivial task.
			
			

		cout<<"Finished step "<<time<<"/"<<number_of_steps<<"\n";
	}



	particles[0].writeDepositToCSV();
	particles[1].writeDepositToCSV();

	cout<<"Simulation finished, writing to csv..."<<"\n";
    //////////////////////////////////////////////////////////////////

	///////////////////COMPUTE ELAPSED TIME///////////////////////////
   	auto endTimer = std::chrono::high_resolution_clock::now();
   	std::chrono::duration<double> duration = endTimer - startTimer;
	double elapsed_seconds = duration.count();
  	std::cout << "Program completed after: " << elapsed_seconds << " seconds" << "\n";
	//////////////////////////////////////////////////////////////////	
	//////////////////END PROGRAM/////////////////////////////////////

	return 0;
}	
