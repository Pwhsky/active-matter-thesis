#ifndef PARTICLE_H
#define PARTICLE_H

struct Particle {
  
	long double x,y;
	double phi;
	bool isHot;
	
	//Phoretic velocity:
	long double vpx,vpy;
};
	
	std::vector<Particle> initialize_stator(); 	//2 hot 2 cold
	std::vector<Particle> initialize_spinner();	//2 hot 4 cold
	std::vector<Particle> initialize_rotator(); 	//2 hot 3 cold
	std::vector<Particle> initialize_migrator();    //1 hot 1 cold
	

#endif // PARTICLE_H
