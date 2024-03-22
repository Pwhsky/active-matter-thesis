import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle
import matplotlib.animation as animation

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
	x = trajectory[100:-1,0]
	coefficients, covariance = np.polyfit(np.log(x),np.log(np.squeeze(MSD[2][100:-1])),1,cov=True)
	errors                   = np.sqrt(np.diag(covariance))

	#from b = 2*log(v) where v is velocity, we can derive velocity with error bars
	#a and b come from coefficients from numpy polyfit.
	#coefficients[0] = 2
	#print(coefficients)
	#print(errors)

	velocity       = np.exp((coefficients[0])/2)
	velocity_error = np.exp((coefficients[0]+errors[0])/2)-np.exp((coefficients[0]-errors[0])/2)

	print(f"velocity = {round(velocity,5)} ± {round(velocity_error,5)}" )
	return velocity,velocity_error
	
##############################################################################################

loop_limit = 4
bounds    = 5
timesteps = 500
nDeposits = 1

velocities    = np.zeros((loop_limit-1,2))
periodicities = np.zeros((loop_limit-1,2))
#Perform parameter sweep over laser periodicities

for i in range(loop_limit):
	periodicity = 1+i*1e-9
	generateNewData(nDeposits, timesteps, periodicity*1e-9,bounds)
	trajectory       = np.genfromtxt("trajectory.csv",delimiter=',',skip_header=3)
	velocity,error   = data_analysis(trajectory)
	velocities[i,0]  = velocity
	velocities[i,1]  = error
	periodicities[i] = periodicity
os.chdir("..")
os.chdir("figures")

plt.plot(velocities,periodicities)
plt.savefig("test.png")
##################

