import matplotlib.pyplot as plt
import numpy as np
import os
import time 
import subprocess
import sys

os.chdir("src")
faceColor = (1,0.964706,0.901961)

limit = 15e-6
font = 19
pi = 3.14159
particleRadius = 2e-6


def generateNewData(nDeposits, timesteps, periodicity,bounds):
	print(f"Producing {timesteps} steps \n" )
	#Compiles and runs the .CPP files
	subprocess.run(["g++","functions.cpp","particle.cpp","brownian_sim.cpp",
					"-o","sim","-Ofast", "-fopenmp","-funroll-all-loops","-march=native"])
	subprocess.run(["./sim",str(nDeposits),str(bounds), str(periodicity),str(timesteps)])

	
def compute_MSD(x):
	n = len(x)
	MSD = np.zeros((n,1))
	for t in range(n):
		square_sum = 0
		for tau in range(n-t):
			square_sum += (x[tau+t]-x[tau])**2
		MSD[t] = square_sum/(n-t)
	return MSD

def data_analysis(trajectory):

	#Trajectory.csv contains time,x,y,z,velocity

	#Rescale to be in units of micrometers squared
	MSD = np.asarray([compute_MSD(trajectory[:,1])*1e12, 
				      compute_MSD(trajectory[:,2])*1e12,
				  	  compute_MSD(trajectory[:,3])*1e12	])

	#Method 1 to obtain slope of MSD, divide MSD with t^2
	#the slopes are in 3 regions of the plot, t[0:20] t[20:400] etc 
	slope_x = np.asarray( [  np.log(MSD[0][20]/MSD[0][1])    /np.log(20), 
							 np.log(MSD[0][400]/MSD[0][20])  /np.log(380/20),
							 np.log(MSD[0][-500]/MSD[0][400])/np.log((len(MSD[0][400:])-500)/380 )  ])
	slope_y = np.asarray( [  np.log(MSD[1][20]/MSD[1][1])    /np.log(20), 
							 np.log(MSD[1][400]/MSD[1][20])  /np.log(380/20),
							 np.log(MSD[1][-500]/MSD[1][400])/np.log((len(MSD[1][400:])-500)/380 )  ])
	slope_z = np.asarray( [  np.log(MSD[2][20]/MSD[2][1])    /np.log(20), 
							 np.log(MSD[2][400]/MSD[2][20])  /np.log(380/20),
							 np.log(MSD[2][-1]/MSD[2][400])  /np.log(len(MSD[2][400:])/380 )  ])

	#Method 2, perform logarithm polyfit (with error estimate):
	x = trajectory[100:-400,0]
	coefficients, covariance = np.polyfit(np.log(x),np.log(np.squeeze(MSD[2][100:-400])),1,cov=True)
	errors                   = np.sqrt(np.diag(covariance))

	#from b = 2*log(v) where v is velocity, we can derive velocity with error bars
	#a and b come from coefficients from numpy polyfit.
	#coefficients[0] = 2
	#print(coefficients)
	#print(errors)

	velocity       = np.exp((coefficients[0])/2)
	velocity_error = np.exp((coefficients[0]+errors[0])/2)-np.exp((coefficients[0]-errors[0])/2)

	#print(f"velocity = {round(velocity,5)} ± {round(velocity_error,5)}" )
	return velocity,velocity_error
	
##############################################################################################


bounds        = 5
timesteps     = 5000
deposit_number = [400,425,450,475,500,
				 525,550,575,600,
				 625,650,675,700,
				 725,750,775,800,
				 825,850,875,900]
velocities = []
errors = []
#Perform parameter sweep over laser periodicities
tic = time.time()
#for i in range(len(deposits)):
#	nDeposits = deposits[i]
#	generateNewData(nDeposits, timesteps, 80000*1e-9,bounds)
#	trajectory       = np.genfromtxt("trajectory.csv",delimiter=',',skip_header=3)
#	velocity,error   = data_analysis(trajectory)
	
#	velocities.append(velocity)
#	errors.append(error)
print(velocities)
print(errors)

deposits = np.asarray(deposit_number)

r = 30e-9
R = 2e-6
particle_volume = (4*pi*R**3)/3
deposit_volume  = (4*pi*r**3)/3

density_FeO   = 5250  #iron oxide density   in Kg/m³
density_Ag    = 19300 #gold density        in kg/m³
density_SiO   = 2320  #silica oxide density in kg/m³

total_iron   = deposits*deposit_volume*density_FeO
total_gold   = deposits*deposit_volume*density_Ag
total_silica = (particle_volume-(deposits*deposit_volume))*density_SiO

mass_ratios  = total_iron/(total_silica+total_iron) 


toc = time.time()
os.chdir("..")
os.chdir("figures")
fig, ax = plt.subplots(figsize=(10, 7))
ax.errorbar(mass_ratios,velocities,errors,capsize=4)

ax.set_ylabel("Velocity",fontsize=20)
ax.set_xlabel(r"$\frac{m_{\rm FeO}}{m_{\rm total}}$ ",fontsize=20)
ax.set_title("Velocity mass relationship",fontsize=20)
plt.savefig("vel_mass.png")


print("Finished after " + str(round(toc-tic)) + " s")
##################

