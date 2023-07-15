#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <random>
#include <omp.h>
#include <fstream>
using namespace std;

const double pi			    = 3.14159265358979323846;
const double temperature            = 300;
const double particleRadius         = 1*pow(10,-6); //1 micrometer sized particle
const double cutoff_distance	    = 8*pow(10,-6);
const double delta_t		    = 0.01;
const double kb      		    = 1.38 * pow(10,-23);
const double viscosity              = 1*pow(10,-3);
const double D_T 	            = kb *temperature / (6.0*pi*viscosity*particleRadius);
const double D_R	   	    = 1;
const double v			    = 2*pow(10,-6);  
const double v_0		    = 10*pow(10,-6); //Phoretic strength
const double box_size		    = pow(10,-4);
const double phoreticTerm	    = v_0 * pow(2*particleRadius,2);
const double trans 		    = sqrt(2.0*D_T*delta_t);
const double rotate 		    = sqrt(2.0*D_R*delta_t);


const double spawnArea		    = box_size/20;




//Random number generator:
	std::random_device rd;
	std::mt19937 gen(rd());

struct Particle{


	//Coords
	double x;
	double y;
	
	//Phoretic velocity directions;
	bool isHot;
	double vpx; 
	double vpy;
	double phi;
}; 


std::vector<Particle> initialize_particles(int nHot, int nCold);
void update_position(Particle &particle);

std::vector<double> phoretic_force(Particle &particle1, Particle &particle2);
void hard_sphere_correction(std::vector<Particle>& particles,int nParticles);

vector<double> getDirection(Particle particle1, Particle particle2);
double getNorm(Particle particle1, Particle particle2);

void boundary_box_correction(std::vector<Particle>& particles,int nParticles);
void writeToCSV(std::vector<std::vector<Particle>> particlesOverTime,int nParticles,int timeSteps);

int main(int argc, char** argv) {
	
	//Take user input: 
	if (argc < 4) {
		std::cout << "Usage: ./sim <nHot> <nCold> <timeSteps>" << std::endl;
		std::cout << "example: \n" << "./sim 2 2 100"<< std::endl;
		return 1;
	}
	
	const int nHot = std::stoi(argv[1]);
	const int nCold = std::stoi(argv[2]);
	const int timeSteps = std::stoi(argv[3]);
	const int nParticles = nHot + nCold;
	
	
	//////////////////
	
	
	//Start a timer
	auto start = std::chrono::high_resolution_clock::now();

	vector<Particle> particles = initialize_particles(nHot,nCold);
	vector<vector<Particle>> particlesOverTime;

   	
	for (int time = 0; time < timeSteps; time++){
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		
		for (int i = 0; i< nParticles; i++){
			double vpx = 0.0;
			double vpy = 0.0;
		
			for(int j = 0; j< nParticles; j++){
				if (i==j) {continue;}
				//Only attracted to hot particles, and not to itself.
				if (particles[j].isHot == true){
					vector<double> totalPhoreticForce = phoretic_force(particles[i],particles[j]);
					vpx += totalPhoreticForce[0];
					vpy += totalPhoreticForce[1];		
				}
			}
			
			particles[i].vpx = vpx;
			particles[i].vpy = vpy;
			update_position(particles[i]);
			
		}
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		boundary_box_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);
		hard_sphere_correction(particles,nParticles);

		particlesOverTime.push_back(particles);
	
		
	}
	writeToCSV(particlesOverTime,nParticles,timeSteps);

   	auto end = std::chrono::high_resolution_clock::now();
   	std::chrono::duration<double> duration = end - start;
   	double elapsed_seconds = duration.count();
  	std::cout << "Simulation completed after: " << elapsed_seconds << " seconds" << std::endl;

	return 0;
}




std::vector<Particle> initialize_particles(int nHot, int nCold) {

	vector<Particle> particles;
	uniform_real_distribution<double> dis(0.0,1.0);
	
	//Assign values to the particles
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
		newParticle.phi = dis(gen)*2.0*pi;
		particles.push_back(newParticle);	
		
	}	
	
	
	
	return particles;
}


void update_position(Particle &particle) {

	normal_distribution<double> normdis(0, 1);
	double W_phi  = normdis(gen);
	double W_x    = normdis(gen);
	double W_y    = normdis(gen);
	double newPhi = particle.phi + W_phi*rotate;
	
	particle.x = particle.x +  v*cos(newPhi) * delta_t + trans*W_x + particle.vpx*delta_t;
	particle.y = particle.y +  v*sin(newPhi) * delta_t + trans*W_y + particle.vpy*delta_t;

}


void hard_sphere_correction(std::vector<Particle> &particles,int nParticles) {


	for (int i = 0 ; i<nParticles-1; i++){
		for(int j = i+1; j<nParticles; j++) {
		
			double centerToCenterDistance =  getNorm(particles[i],particles[j]);
			double overlap = 2.0*particleRadius - centerToCenterDistance;
			
			if (overlap > 0 && overlap != 2.0*particleRadius){
			
				double distanceToMove = overlap/2.0;
				vector<double> direction = getDirection(particles[i],particles[j]);
				
				particles[i].x = particles[i].x + distanceToMove*direction[0];
				particles[i].y = particles[i].y + distanceToMove*direction[1];	
				
				particles[j].x = particles[j].x - distanceToMove*direction[0];
				particles[j].y = particles[j].y - distanceToMove*direction[1];	
				
			
			}						
		}
	}
}



//Compute the velocity of the phoretic interaction for a given particle
std::vector<double> phoretic_force(Particle &particle1, Particle &particle2){
	double vpx = 0.0;
	double vpy = 0.0;
	
	double particleDistance = getNorm(particle1,particle2);
			
	if (particleDistance < cutoff_distance && particleDistance !=0) {
		double directionx =  (particle2.x-particle1.x)/particleDistance;
		double directiony =  (particle2.y-particle1.y)/particleDistance;
		vpx = (phoreticTerm/(pow(particleDistance,2)))*directionx;
		vpy = (phoreticTerm/(pow(particleDistance,2)))*directiony;
		
		return {vpx,vpy};
	} else {
		return {vpx,vpy};
	}

}


vector<double> getDirection(Particle particle1, Particle particle2){
	//Return normalized direction, in x and y coordinates.
	double norm = sqrt(pow(particle1.x-particle2.x,2) + pow(particle1.y-particle2.y,2));
	double x = (particle1.x-particle2.x)/norm;
	double y = (particle1.y-particle2.y)/norm;
	
	vector<double> direction = {x,y};

	
	return direction;
}


double getNorm(Particle particle1, Particle particle2){
	double norm = sqrt(pow(particle1.x-particle2.x,2) + pow(particle1.y-particle2.y,2));
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
  		else
  		{
  			coldFile << particle.x << "," << particle.y << std::endl;
  		}
    	}	   
    
    }
    	hotFile.close();
   	coldFile.close();
}




