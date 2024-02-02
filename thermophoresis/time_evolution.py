import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import concurrent.futures
import os
import time 
import subprocess
import sys

os.chdir("src")
from matplotlib.patches import Circle
from cython_functions import histogram2d_cython, gradient_cython

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
	generateData    	       = str(sys.argv[2])
	resolution 	        	   = 300 
	nDeposits                  = str(sys.argv[1])
	return resolution,nDeposits,generateData
	
def generateNewData():
	#Compiles and runs the .CPP file
	print("Generating new data...\n")
	subprocess.run(["g++","functions.cpp","particle.cpp","brownian_sim.cpp","-o","sim","-Ofast", "-fopenmp" , "-funroll-all-loops"])
	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4]])

def loadData():
	df = pd.read_csv("gradient.csv",engine="pyarrow")
	df.sort_values(by=['z'])

	with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    		#Parallel data loading
    		futures = [
        	executor.submit(pandasToNumpy, df['x']),
        	executor.submit(pandasToNumpy, df['y']),
      		executor.submit(pandasToNumpy, df['z']),

       		]
	x    = futures[0].result()
	y    = futures[1].result()
	z    = futures[2].result()

	return x,y,z
	
def generateFigure(imageBounds):
	
	depositDF = pd.read_csv('deposits.csv')
	initial   = pd.read_csv("initpositions.csv",engine="pyarrow")
	final     = pd.read_csv("finalpositions.csv",engine="pyarrow")

	fig, ax     = plt.subplots(1, 1, figsize=(13, 7))
	ax.set_title("Final orientation of particle")


	xlim = [final['x'][0] - 2e-5, final['x'][0] + 2e-5]
	ylim = [final['z'][0] - 2e-5, final['z'][0] + 2e-5]

	circle1 = Circle((final['x'][0], final['z'][0]), 2e-6)
	circle1.set(fill=False, alpha=0.5)

	circle2 = Circle((final['x'][1], final['z'][1]), 2e-6)
	circle2.set(fill=False, alpha=0.5)

	ax.add_patch(circle1)
	ax.add_patch(circle2)
	ax.scatter(depositDF['x'][:],depositDF['z'][:],s=10)
	#ax.scatter(lastPos[0],lastPos[1],label="Particle Center")
	for i in range(len(initial['x'])):
		ax.scatter(final['x'][i],final['z'][i]				,label="Starting position")
		ax.plot(initial['x'][i],final['z'][i],linestyle="--",label="Trajectory")	
	
	
	#ax.hlines(lastPos[1],lastPos[0],lastPos[0]+2e-6)
	ax.axis("equal")
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	ax.set_xlabel('X (m)',fontsize=15)
	ax.set_ylabel('Z (m)',fontsize=15)
	ax.grid(True)
	ax.legend()



resolution,nDeposits,generateData = parseArgs()
if (generateData == "true"):
	generateNewData()
else:
	print("No new data will be generated, starting plotting... \n")
	
tic = time.time()
print("Starting plotting...")
x,y,z= loadData()

#Saves the figure in the figures directory, then goes back to generate a new figure and 
#then saves in the figures directory once more.
generateFigure(imageBounds)
os.chdir("..")
os.chdir("figures")

plt.savefig("time_evolution.png")

toc = time.time()
#print("Plotting finished after " + str(round(toc-tic)) + " s")
