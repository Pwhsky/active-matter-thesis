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
deposits 	  = [200,250,275,300,375,
				 400,425,450,475,500,550,575,
				 600,625,650,675,700]
velocities    = [2.331652584771042, 2.3921237210005333, 2.508021240611387, 2.4391883893443795,   2.4830034948199877, 
				2.4458827912838013, 2.4700890866227843, 2.636457814022954, 2.526986576475101, 2.535214644356422, 2.581321191753488, 2.642019075707256, 2.5749741316427315, 2.583056302884516, 2.6214465830135714, 2.6277173448088655, 2.65464558807304]
errors  	  = [0.006084998729259716, 0.005993665345243482, 0.003646279993120949, 0.0032192750955553073, 0.003999915647349361, 
				0.00556534189141944, 0.002807099718984407, 0.0021604297117203686, 0.003160881020285, 0.0030998901066414675, 0.003121191055168726, 0.0022468626906535505, 0.0031772872563911037, 0.0014760774007047672, 0.002336761070151905, 0.001990607739728034, 0.0012425598387206804]
#errors 		  = [0.2,0.4,0.8,1.6, 3.2, 6.4,12.8, 25,6.4, 51.2,102.4]
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

deposits = np.asarray(deposits)
r = 30e-9
R = 2e-6
particle_volume = (4*pi*R**3)/3
deposit_volume  = (4*pi*r**3)/3
ratios = deposits*(deposit_volume/particle_volume)


toc = time.time()
os.chdir("..")
os.chdir("figures")
fig, ax = plt.subplots(figsize=(10, 7))
ax.errorbar(ratios,velocities,errors,capsize=4)

ax.set_ylabel("Velocity",fontsize=20)
ax.set_xlabel(r"$\frac{V_{iron}}{V_{particle}}$ ",fontsize=20)
ax.set_title("Velocity vs volume ratio",fontsize=20)
plt.savefig("test.png")


print("Finished after " + str(round(toc-tic)) + " s")
##################

