#include <cmath>
#include <vector>
#include "particle.h"
using namespace std;



	const double particleRadius         = 1*pow(10,-6);
	const double box_size		    = 100*pow(10,-6);
	const double spawnArea		    = box_size/45;

//This file contains some functions to place pre-defined cluster configurations at the center of the arena, 
//some functions are more easy to read than others...


std::vector<Particle> initialize_stator() {

	vector<Particle> particles;
	
	//Assign values to the particles and then add to the particles list.
	//First hot particles
	for (int i = 1; i < 3; i++) {
		Particle newParticle;
		newParticle.isHot = true;

		newParticle.x =  (2*i*particleRadius) + box_size/2- 3*particleRadius;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
	}	
	
	for (int i = 1; i < 3; i++) {
		Particle newParticle;
		newParticle.isHot = false;

		newParticle.x =  box_size/2;		
		newParticle.y =  (2*i*particleRadius) + box_size/2 - 3*particleRadius;
		newParticle.vpx,newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);	
		
	}	
	return particles;
}


//place one spinner, 2 hot 4 cold
std::vector<Particle> initialize_spinner() {

	vector<Particle> particles;

		Particle newParticle;
		
		newParticle.isHot = false;
		newParticle.x =  box_size/2 -2*particleRadius;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
		
		
		newParticle.isHot = true;
		newParticle.x =  box_size/2;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
		
		newParticle.isHot = false;
		newParticle.x =  box_size/2 +2*particleRadius;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);

		newParticle.isHot = true;
		newParticle.x =  box_size/2 +particleRadius;		
		newParticle.y =  box_size/2 - particleRadius;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
		
		newParticle.isHot = false;
		newParticle.x =  box_size/2 -particleRadius;		
		newParticle.y =  box_size/2 - particleRadius;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
		
		newParticle.isHot = false;
		newParticle.x =  box_size/2 +3*particleRadius;		
		newParticle.y =  box_size/2 - particleRadius;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);


	return particles;
}

//place one rotator,   2 hot 3 cold
std::vector<Particle> initialize_rotator() {

	vector<Particle> particles;

		Particle newParticle;
		
		newParticle.isHot = false;
		newParticle.x =  box_size/2 -2*particleRadius;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
		
		
		newParticle.isHot = true;
		newParticle.x =  box_size/2;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
		
		newParticle.isHot = true;
		newParticle.x =  box_size/2 +1.8*particleRadius;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);

		
		newParticle.isHot = false;
		newParticle.x =  box_size/2 -particleRadius;		
		newParticle.y =  box_size/2 - particleRadius;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
		
		newParticle.isHot = false;
		newParticle.x =  box_size/2 +particleRadius;		
		newParticle.y =  box_size/2 - particleRadius;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);


	return particles;
}

//place one rotator,   2 hot 3 cold
std::vector<Particle> initialize_migrator() {

	vector<Particle> particles;
	
		Particle newParticle;
		
		newParticle.isHot = false;
		newParticle.x =  box_size/2;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
		
		newParticle.isHot = true;
		newParticle.x =  box_size/2;		
		newParticle.y =  box_size/2;
		newParticle.vpx, newParticle.vpy, newParticle.phi = 0.0;
		particles.push_back(newParticle);
	


	return particles;
}

