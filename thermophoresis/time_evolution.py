import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle,Wedge

os.chdir("src")



##IMPORTANT:
#The units on the colorbar are kelvin divided by the step-size of the finite differences.
#When resolution = 200 (for most simulations) that equals 25 nanometers!
#Convert accordingly.

font = 19
pi = 3.14159
particleRadius = 2e-6

imageBounds 	   = float(sys.argv[3])*1e-6
spatialPeriodicity = float(sys.argv[4])*1e-9
periodicity 	   = float(sys.argv[4])/1000

def pandasToNumpy(column):
    return column.to_numpy()

def parseArgs():
	resolution 	        	   = 300 
	nDeposits                  = str(sys.argv[1])
	return resolution,nDeposits
	
def generateNewData():
	#Compiles and runs the .CPP file
	print(f"Producing {sys.argv[2]} steps \n" )
	#Forbidden compiler flags below, if it doesn't break anything then leave them in (10% speed boost)
	subprocess.run(["g++","functions.cpp","particle.cpp","brownian_sim.cpp","-o","sim","-Ofast", "-fopenmp", "-fno-strict-overflow","-fno-unsafe-loop-optimizations", "-fno-strict-aliasing","-funroll-all-loops"])
	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4],sys.argv[2]])


	
def generateFigure(imageBounds):
	
	#depositDF = pd.read_csv('deposits.csv')
	particle_1     = pd.read_csv("particle_1.csv",engine="pyarrow")
	particle_2     = pd.read_csv("particle_2.csv",engine="pyarrow")

	fig, ax     = plt.subplots(1, 1, figsize=(13, 7))
	ax.set_title(f"Time evolution after {sys.argv[2]} timesteps.")
	
	xlim = [particle_1.iloc[-1]['x'] - 3e-5, particle_1.iloc[-1]['x'] + 3e-5]
	ylim = [particle_1.iloc[-1]['x'] - 3e-5, particle_1.iloc[-1]['x'] + 3e-5]

	circle1 = Circle((particle_1.iloc[-1]['x'], particle_1.iloc[-1]['z']), 2e-6)
	wedge1  =  Wedge((particle_1.iloc[-1]['x'], particle_1.iloc[-1]['z']), 2e-6, 0, 180, color='red', alpha=0.6)
	circle1.set(fill=False, alpha=0.5)

	circle2 = Circle((particle_2.iloc[-1]['x'], particle_2.iloc[-1]['z']), 2e-6)
	wedge2  =  Wedge((particle_2.iloc[-1]['x'], particle_2.iloc[-1]['z']), 2e-6, 0, 180, color='red', alpha=0.6)
	circle2.set(fill=False, alpha=0.5)

	ax.add_patch(circle1)
	ax.add_patch(wedge1)
	ax.add_patch(circle2)
	ax.add_patch(wedge2)
	#ax.scatter(depositDF['x'][:],depositDF['z'][:],s=10)
	#ax.scatter(lastPos[0],lastPos[1],label="Particle Center")


	ax.plot(particle_1['x'][:],particle_1["z"][:],label="Trajectory 1",color='black')
	ax.plot(particle_2['x'][:],particle_2["z"][:],label="Trajectory 2",color='blue')
	
	#ax.hlines(lastPos[1],lastPos[0],lastPos[0]+2e-6)
	ax.axis("equal")
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	ax.set_xlabel('X (m)',fontsize=15)
	ax.set_ylabel('Z (m)',fontsize=15)
	ax.legend()
	os.chdir("..")
	os.chdir("figures")
	plt.savefig("time_evolution.png")
	

def plotDeposits():

	os.chdir("..")
	os.chdir("src")
	df          = pd.read_csv('deposits.csv')
	fig, ax     = plt.subplots(1, 1, figsize=(13, 7))
	ax = fig.add_subplot(projection='3d')
	ax.view_init(azim=50, elev=10)
	ax.set_box_aspect([1, 1, 1])
	#ax.scatter3D(df['x'],df['y'],df['z'])
	os.chdir("..")
	os.chdir("figures")
	plt.savefig("deposits.png")


resolution,nDeposits = parseArgs()

generateNewData()

tic = time.time()
print("Starting plotting...")

#Saves the figure in the figures directory, then goes back to generate a new figure and 
#then saves in the figures directory once more.
generateFigure(imageBounds)
#plotDeposits()
toc = time.time()
#print("Plotting finished after " + str(round(toc-tic)) + " s")
