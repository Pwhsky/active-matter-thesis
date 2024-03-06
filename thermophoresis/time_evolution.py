import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle

os.chdir("src")



##IMPORTANT:
#The units on the colorbar are kelvin divided by the step-size of the finite differences.
#When resolution = 200 (for most simulations) that equals 25 nanometers!
#Convert accordingly.

font = 19
pi = 3.14159
particleRadius = 2e-6

nDeposits          = str(sys.argv[1])
imageBounds 	   = float(sys.argv[3])*1e-6
spatialPeriodicity = float(sys.argv[4])*1e-9

def generateNewData():
	#Compiles and runs the .CPP files
	print(f"Producing {sys.argv[2]} steps \n" )
	subprocess.run(["g++","functions.cpp","particle.cpp","brownian_sim.cpp","-o","sim","-Ofast", "-fopenmp","-funroll-all-loops"])
	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4],sys.argv[2]])

def generateFigure(imageBounds):

	particle_1     = pd.read_csv("particle_1.csv",engine="pyarrow")
	#particle_2     = pd.read_csv("particle_2.csv",engine="pyarrow")
	deposits       = pd.read_csv('deposits.csv')

	fig, ax     = plt.subplots(1, 1, figsize=(10, 7))
	ax.set_title(f"Time evolution after {sys.argv[2]} timesteps.")
	
	
	limit = 25e-6
	xlim = [- limit,  limit]
	ylim = [- 5*limit,  limit]

	#Create circles for particles:
	circle1 = Circle((particle_1.iloc[-1]['x'], particle_1.iloc[-1]['z']), 2e-6)
	#circle2 = Circle((particle_2.iloc[-1]['x'], particle_2.iloc[-1]['z']), 2e-6)
	circle1.set(fill=False, alpha=0.5)
	#circle2.set(fill=False, alpha=0.5)
	##############################

	#Generate contour plot of laser profile:
	x   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
	y   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
	X,Y = np.meshgrid(x,y)	
	Z = (1+ (np.cos(2*pi*X/(spatialPeriodicity))))
	ax.contourf(X,Y,Z,400,alpha=1,zorder=-1)
	######################################

	#Plot, place, and draw the deposits, circles, trajectory:
	ax.plot(particle_1['x'][:],particle_1["z"][:],label="Trajectory 1",color='black',zorder=0)
	#ax.plot(particle_2['x'][:],particle_2["z"][:],label="Trajectory 2",color='blue',zorder=0)
	ax.scatter(deposits['x'],deposits['z'],color='red',s=10,alpha=1,zorder=0)
	ax.add_patch(circle1)
	#ax.add_patch(circle2)

	###########################################################

	ax.axis("equal")
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	ax.set_xlabel('X (m)',fontsize=15)
	ax.set_ylabel('Z (m)',fontsize=15)
	ax.legend()

	#Save figure in figures directory:
	os.chdir("..")
	os.chdir("figures")
	plt.savefig("time_evolution.png")
	#################################

def plotDeposits():

	os.chdir("..")
	os.chdir("src")
	deposits    = pd.read_csv('deposits.csv')
	fig, ax     = plt.subplots(1, 1, figsize=(13, 7))
	ax.axis("equal")
	ax.scatter(deposits['x'],deposits['z'])
	os.chdir("..")
	os.chdir("figures")
	plt.savefig("deposits.png")



generateNewData()
print("Starting plotting...")
generateFigure(imageBounds)

#plotDeposits()