#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <random>
#include <omp.h>
using namespace std;

const double pi			    = 3.14159265358979323846;
const double temperature            = 300;
const double viscosity              = 1*pow(10,-3);
const double particleRadius         = 1 * pow(10,-6); //1 micrometer sized particle
const double kb      		    = 1.38 * pow(10,-23);
const double D_T 	            = kb *temperature / (6*pi*viscosity);
const double D_R	   	    = 1;
const double v			    = 2*pow(10,-6);
const double v_0		    = 50*pow(10,-6);
const double delta_t		    = 0.03;
const double box_size		    = 100*pow(10,-6);
const double cutoff_distance	    = 8*particleRadius;
const double phoreticTerm	    = v_0 * pow(particleRadius,2);
const double trans 		    = sqrt(2.0*D_T*delta_t);
const double rotate 		    = sqrt(2.0*D_R*delta_t);
const int    nParticles		    = 1000;
const int    timeSteps		    = 100;

//Random number generator:
	std::random_device rd;
	std::mt19937 gen(rd());

struct Particle{

	double x;
	double y;
	
	//Phoretic velocities;
	double vpx; 
	double vpy;
	
	double phi;
}; 

std::vector<double> phoretic_force(Particle &particle1, Particle &particle2);
std::vector<Particle> initialize_particles(int nParticles);
void update_position(Particle &particle);


int main() {
	auto start = std::chrono::high_resolution_clock::now();

	vector<Particle> particles = initialize_particles(nParticles);
	
	
	for (int time = 0; time < timeSteps; time++){
		
	
		for (int i = 0; i< nParticles; i++){
			double vpx = 0;
			double vpy = 0;
		
			for(int j = 0; j< nParticles; j++){
				vector<double> totalPhoreticForce = phoretic_force(particles[i],particles[j]);
				vpx += totalPhoreticForce[0];
				vpy += totalPhoreticForce[1];		
			}

			particles[i].vpx = vpx;
			particles[i].vpy = vpy;
			update_position(particles[i]);
		}
		//Check for overlaps and perform hard sphere correction where needed:
		
			
		
		
	}
	
	
	
	
	
	
	
	
	
    auto end = std::chrono::high_resolution_clock::now();

   
    std::chrono::duration<double> duration = end - start;
    double elapsed_seconds = duration.count();

    
    std::cout << "Elapsed time: " << elapsed_seconds << " seconds" << std::endl;

	

	return 0;
}




std::vector<Particle> initialize_particles(int nParticles) {
	vector<Particle> particles;
	
	
	uniform_real_distribution<double> dis(0.0,1.0);
	
	//Assign values to the particles
	for (int i = 0; i<nParticles; i++) {
		Particle newParticle;
		newParticle.x = dis(gen)*box_size;
		newParticle.y = dis(gen)*box_size;
		newParticle.vpx = 0.0;
		newParticle.vpy = 0.0;
		newParticle.phi = dis(gen)*2*pi;
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
	particle.y = particle.y +  v*cos(newPhi) * delta_t + trans*W_y + particle.vpy*delta_t;

}
void hard_sphere_correction(std::vector<Particle> particles) {
	for (int i = 0 ; i<nParticles; i++){
		for(int j = i; j<nParticles; j++) {
			
			
			
		}
	}

}

std::vector<double> phoretic_force(Particle &particle1, Particle &particle2){
	double vpx = 0.0;
	double vpy = 0.0;
	
	double particleDistance = sqrt( pow(particle2.x-particle1.x,2) + pow(particle2.y-particle1.y,2) );
			
	if (particleDistance < cutoff_distance && particleDistance !=0) {
		double directionx =  (particle2.x-particle1.x)/particleDistance;
		double directiony =  (particle2.y-particle1.y)/particleDistance;
		vpx = phoreticTerm/(pow(particleDistance,2))*directionx;
		vpy = phoreticTerm/(pow(particleDistance,2))*directiony;
		
		return {vpx,vpy};
	} else {
	
		return {vpx,vpy};
	}


}




/*class Particle {
public:
	Particle(double x, double y)
	: position({x,y}) {}
	
	double getX() const {return position.x;}
	double getY() const {return position.y;}
	
	void setX(float x) { position.x = x; }
        void setY(float y) { position.y = y; }
	
private:
    struct {
    	double x;
    	double y;
    } position;
}; */
