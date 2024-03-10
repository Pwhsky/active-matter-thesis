/*
Alex Lech 2023	

This source code contains the time evolution for simulating thermophoresis.
*/

#include <chrono>
#include <fstream>
#include <cmath>
#include <vector>
#include "particle.h"
double compute_total_velocity(Particle &particle);

using namespace std;
	int     number_of_steps,nDeposits;
 	double  bounds, lambda, stepSize, dv, dl; 
	bool    onlyTangential = false;
	vector<Point> globalDeposits;


int main(int argc, char** argv) {
	   
	auto startTimer = std::chrono::high_resolution_clock::now();
	bounds    		= 100  * pow(10,-6); 
	stepSize  		= bounds/(double)(300);		      //Step size, based off of bounds parameter Set to 100 for now
	nDeposits 		= stof(argv[1]);				  //number of deposits to initialize
	lambda	 	    = stold(argv[3])  * pow(10,-9);   //Spatial periodicity
	number_of_steps = (int)stof(argv[4]);
    dv	      		= stepSize*stepSize*stepSize;     //volume element for integral


	std::ofstream write_deposit_data("deposits.csv");
	std::ofstream write_trajectory_data("trajectory.csv");
	//std::ofstream p2("particle_2.csv");

	write_deposit_data      << "x,y,z"<<"\n";
	write_trajectory_data   << "t,x,y,z,v"<<"\n";
	//p2   << "x,y,z"<<"\n";

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

	
	//Write initial positions to file:
	write_trajectory_data << 0.0 <<"," 
						  << particles[0].center.x << ","
	   					  << particles[0].center.y << "," 
	   					  << particles[0].center.z << ","
						  <<                   0.0 << "\n";
	for(int i = 0; i<globalDeposits.size(); i++){
		write_deposit_data << globalDeposits[i].x << "," 
		   				   << globalDeposits[i].y << "," 
		   			 	   << globalDeposits[i].z <<"\n";
	}
	
//	p2 << particles[1].center.x << "," << particles[1].center.y << "," << particles[1].center.z << "\n";

	for(int time = 0; time < number_of_steps; time ++){ 


		for(auto &particle:particles){
			particle.getKinematics(linspace,thickness,dl,globalDeposits,lambda,dv);
		}

		for(auto &particle:particles){
			particle.rotation_transform();
			particle.update_position();
			particle.brownian_noise();
		}
		hard_sphere_correction(particles);


		//Write deposits to csv file:
		for(int i = 0; i<globalDeposits.size(); i++){
			write_deposit_data << globalDeposits[i].x << "," 
			   				   << globalDeposits[i].y << "," 
			   				   << globalDeposits[i].z <<"\n";
		}
		globalDeposits = update_globalDeposits(particles);
		
		//Write velocity and displacements to csv file
		double total_vel = compute_total_velocity(particles[0]);
		write_trajectory_data << time*0.01 << ","
		                      << particles[0].center.x << ","
							  << particles[0].center.y << ","
							  << particles[0].center.z << "," 
							  << total_vel             << "\n";
		//p2 << particles[1].center.x << "," << particles[1].center.y << "," << particles[1].center.z << "\n";
			
		cout<<"Finished step "<<time<<"/"<<number_of_steps<<"\n";
	}

	write_trajectory_data.close();
	write_deposit_data.close();

	cout<<"Simulation finished, writing to csv..."<<"\n";
    //////////////////////////////////////////////////////////////////

	///////////////////COMPUTE ELAPSED TIME///////////////////////////
   	auto endTimer = std::chrono::high_resolution_clock::now();
   	std::chrono::duration<double> duration = endTimer - startTimer;
	double elapsed_seconds = duration.count();
  	std::cout << "Simulation finished after: " << elapsed_seconds << " seconds" << "\n";
	//////////////////////////////////////////////////////////////////	
	//////////////////END PROGRAM/////////////////////////////////////

	return 0;
}	

double compute_total_velocity(Particle &particle){
	double total_vel = 0;
		for(int i = 0; i<3; i++){
			total_vel += particle.velocity[i]*particle.velocity[i];
		}
		total_vel = sqrt(total_vel);
	return total_vel;
}