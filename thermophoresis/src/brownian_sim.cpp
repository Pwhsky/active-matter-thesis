/*
Alex Lech 2023	

This source code contains the computation for the temperature gradient.
compute_temperature.cpp contains the computation for the temperature increase.
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
			particle.updatePosition();
			
			hard_sphere_correction(particles);


			p1 << particles[0].center.x << "," << particles[0].center.y << "," << particles[0].center.z << "\n";
			p2 << particles[1].center.x << "," << particles[1].center.y << "," << particles[1].center.z << "\n";
			//TODO: make exporting particle positions a trivial task.
			
		}
			
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

/*
double central_difference(double x_back,double x_forward,double y_back,double y_forward, double z_back, double z_forward, vector<Point> deposits){
	double back   		= integral(x_back,y_back,z_back,deposits);
	double forward		= integral(x_forward,y_forward,z_forward,deposits);
	return (forward - back)/(2*dl);
}

inline double integral(double _x,double _y,double _z,std::vector<Point> deposits){
	//absorbtionTerm will compute the absorbed ammount of power from the laser
	//ContributionSum will sum up contributions from all deposits
	//Finally, the contributionSum is scaled with volume element dv and divided with constants												
	double laserPower	           = I0 + I0*cos(twoPi*(_x)/lambda);	
	double absorbtionTerm          = laserPower*depositArea/(volumePerDeposit);
	double contributionSum 		   = 0.0;
	
	//Since the values scale with the inverse square distance.
    	for (size_t i = 0; i < deposits.size(); i++){

    		double inverse_squareroot_distance = 1.0/sqrt(pow(_x-deposits[i].x,2)+
														  pow(_y-deposits[i].y,2)+
														  pow(_z-deposits[i].z,2));
			contributionSum +=  inverse_squareroot_distance;
		}
    return contributionSum*absorbtionTerm*dv/(4*pi*waterConductivity); 
}
*/
