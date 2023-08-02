#include <iostream>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <fstream>
#include "particle.h"

using namespace std;
std::random_device rd;
std::mt19937 gen(rd());


//Constants:
	const double pi			    = 3.14159265358979323846;
	const double temperature            = 300;
	const double particleRadius         = 1*pow(10,-6); 			   //1 micrometer sized particle
	const double cutoff_distance	    = 8*pow(10,-6);
	const double delta_t		    = 0.005;
	const double kb      		    = 1.380649 * pow(10,-23);
	const double viscosity              = 1*pow(10,-3);
	const double stokesCoefficient	    = 6.0*pi*viscosity*particleRadius;
	const double D_T 	            = kb *temperature / stokesCoefficient; //Translational diffusion
	const double v_0		    = 10*pow(10,-6);			   //Phoretic strength
	const double box_size		    = 100*pow(10,-6);
	const long double phoreticTerm	    = 4*v_0 * pow(particleRadius,2);
	const long double trans 	    = sqrt(2.0*D_T*delta_t)*0.1;
	const double v			    = particleRadius*10;
	const double spawnArea		    = box_size/45;



//Random initialization
	std::vector<Particle> initialize_particles(int nHot, int nCold);



void update_position	     (Particle &particle,std::vector<Particle> &particles, int nParticles);
void phoretic_force          (std::vector<Particle>& particles,int nParticles);
void hard_sphere_correction  (std::vector<Particle>& particles,int nParticles);
void boundary_box_correction (std::vector<Particle>& particles,int nParticles);
void writeToCSV		     (std::vector<std::vector<Particle>> particlesOverTime,int nParticles,int timeSteps);

vector<long double> getDirection(Particle &particle1, Particle &particle2);
double 		    getNorm(Particle &particle1, Particle &particle2);




int main(int argc, char** argv) {

	auto startTimer = std::chrono::high_resolution_clock::now();
	
	/////////////////////////Take user input: /////////////////
		if (argc < 3) {
			std::cout << "Usage: ./sim <clusterType>  <timeSteps>" << std::endl;
			std::cout << "example: \n" << "./sim spinner 100"<< std::endl;
			return 1;
		}
		std::string clusterType = argv[1];
		const int    timeSteps   = std::stoi(argv[2]);
				
	///////////////////Initialize the cluster///////////////
		vector<vector<Particle>> particlesOverTime;
		vector<Particle> particles; 
		if       (clusterType == "spinner"){
			particles = initialize_spinner();
		}else if (clusterType == "rotator") {
		        particles = initialize_rotator();
		}else if (clusterType == "stator") {
			particles = initialize_stator();
		}else if (clusterType == "migrator"){
			particles = initialize_migrator();
		}
		const int nParticles = particles.size();
	/////////////////////////////////////////////////////
	
	
   	
   	//simulate over time:
	for (int time = 0; time < timeSteps; time++){
	
		////////////////////Phoretic interaction///////////////
			phoretic_force(particles,nParticles);
			for (int i = 0; i<nParticles;i++){
				update_position(particles[i],particles,nParticles);
			}
		///////////////////Hard Sphere Correction////////////
			hard_sphere_correction( particles, nParticles);
			hard_sphere_correction( particles, nParticles);
			hard_sphere_correction( particles, nParticles);
			hard_sphere_correction( particles, nParticles);
		//////////////////////////////////////////////////////
		
		particlesOverTime.push_back(particles);		
	}
	writeToCSV(particlesOverTime,nParticles,timeSteps);


	///////////////////Compute elapsed time/////////////////////////
   		auto endTimer = std::chrono::high_resolution_clock::now();
   		std::chrono::duration<double> duration = endTimer - startTimer;
   		double elapsed_seconds = duration.count();
  		std::cout << "Simulation completed after: " << elapsed_seconds << " seconds" << std::endl;
	////////////////////////////////////////////////////////////////////////////	
	
	return 0;
}




std::vector<Particle> initialize_particles(int nHot, int nCold) {
	//Randomly initializes particles in the center of the arena
	
	vector<Particle> particles;
	uniform_real_distribution<double> dis(0.0,1.0);
	
	for (int i = 0; i< (nHot + nCold); i++) {
		Particle newParticle;
		
		if(i < nHot){
			newParticle.isHot = true;
		}else{
			newParticle.isHot = false;
		}		
		newParticle.x = dis(gen)*spawnArea + box_size/2;
		newParticle.y = dis(gen)*spawnArea + box_size/2;
		newParticle.vpx = 0.0;
		newParticle.vpy = 0.0;
		newParticle.phi = 0.0;
		particles.push_back(newParticle);			
	}	
	return particles;
}


void update_position(Particle &particle,std::vector<Particle> &particles, int nParticles) {
	
	normal_distribution<double> normdis(0.0, 1.0);
	long double W_x    = normdis(gen);
	long double W_y    = normdis(gen);
	double W_phi  = normdis(gen);
	
	if (particle.isHot == true){
		particle.x = particle.x + W_x*trans + particle.vpx*delta_t;
		particle.y = particle.y + W_y*trans + particle.vpy*delta_t;
	} else {
	
		particle.x = particle.x  + particle.vpx*delta_t + W_x*trans;
		particle.y = particle.y  + particle.vpy*delta_t + W_y*trans;
	}
	
}


void hard_sphere_correction(std::vector<Particle> &particles,int nParticles) {

	//Perform synchronous hard sphere correction using a temp.
	std::vector<Particle> tempParticles = particles;
	for (int i = 0; i<nParticles; i++){
		for(int j = 0; j<nParticles; j++) {
			long double centerToCenterDistance =  getNorm(particles[i],particles[j]);
			long double overlap =   2.0*particleRadius- centerToCenterDistance;
			
			if (overlap > 0 && overlap != 0 && i !=j ){
			
				long double distanceToMove = overlap/6.0;
				vector<long double> direction = getDirection(particles[i],particles[j]);				
		
				tempParticles[i].x -=  distanceToMove*direction[0];
				tempParticles[i].y -=  distanceToMove*direction[1];	
				tempParticles[j].x +=  distanceToMove*direction[0];
				tempParticles[j].y +=  distanceToMove*direction[1];	
			}							
		}
	}
	particles = tempParticles;
}



void phoretic_force(std::vector<Particle>& particles,int nParticles){


	//calculate total phoretic force for each particle	
	for (int i = 0; i< nParticles; i++){
		particles[i].vpx = 0.0;
		particles[i].vpy = 0.0;
		
		for(int j = 0; j< nParticles; j++){	
			if (particles[j].isHot == true && i != j){
				double particleDistance = getNorm(particles[i],particles[j]);
				
				if (particleDistance < cutoff_distance && particleDistance != 0) {

					double directionx  =  (particles[i].x-particles[j].x)/particleDistance;
					double directiony  =  (particles[i].y-particles[j].y)/particleDistance;
					particles[i].vpx  -= (phoreticTerm/(pow(particleDistance,2)))*directionx;		
					particles[i].vpy  -= (phoreticTerm/(pow(particleDistance,2)))*directiony;
				}
	
			}		
		}	
	}
}


vector<long double> getDirection(Particle &particle1, Particle &particle2){
	//Return normalized direction, in x and y coordinates.
	
	long double norm = getNorm(particle1,particle2);
	long double x = (particle2.x-particle1.x)/norm;
	long double y = (particle2.y-particle1.y)/norm;
	vector<long double> direction = {x,y};

	return direction;
}


double getNorm(Particle &particle1, Particle &particle2){
	double xdiff = particle1.x-particle2.x;
	double ydiff = particle1.y-particle2.y;

	double norm = sqrt( pow(xdiff,2) + pow(ydiff,2));
	return norm;
}


void boundary_box_correction(vector<Particle> &particles,int nParticles) {

	for (int i = 0; i<nParticles; i++){
		if(particles[i].x > box_size){
       			particles[i].x = particles[i].x - box_size;
       		}     
       		else if (particles[i].x < 0){
        		particles[i].x = particles[i].x + box_size;
       		}
       		if(particles[i].y > box_size){
       			particles[i].y = particles[i].y - box_size;
       		}     
       		else if (particles[i].y < 0){
        		particles[i].y = particles[i].y + box_size;
       		}
	}
}






void writeToCSV(std::vector<std::vector<Particle>> particlesOverTime,int nParticles,int timeSteps) {
	ofstream hotFile("hotParticles.csv");
	ofstream coldFile("coldParticles.csv");
   	hotFile << "x,y" << std::endl;
   	coldFile << "x,y" << std::endl;
   	
    // Write particle positions
	for (int t = 0;t<timeSteps; t++) {
    		vector<Particle> particles = particlesOverTime[t];
    
    		for (const auto& particle : particles) {
    
  			if (particle.isHot == true){
  				hotFile << particle.x << "," << particle.y << std::endl;
  			}
  			else{
  				coldFile << particle.x << "," << particle.y << std::endl;
  			}
    		}	   
    
   	}
    	hotFile.close();
   	coldFile.close();
}

