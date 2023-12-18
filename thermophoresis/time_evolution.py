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
circle1 = Circle((0, 0), 2e-6)
circle1.set(fill=False, linestyle='--', alpha=0.2)
circle2 = Circle((0,0), 2e-6)
circle2.set(fill=False,linestyle='--',alpha=0.2)

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
	subprocess.run(["g++","functions.cpp","particle.cpp","compute_gradient.cpp","-o","sim","-Ofast", "-fopenmp" , "-funroll-all-loops"])
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
	positions = pd.read_csv("positions.csv",engine="pyarrow")
	fig, ax     = plt.subplots(1, 2, figsize=(22, 5))


	ax[0].plot(positions['x'][:],positions['z'][:])

	positions_reset = positions.reset_index(drop=True)

	lastPos = [positions_reset['x'].iloc[-1],positions_reset['z'].iloc[-1]]
	xlim = [lastPos[0] - 2e-6, lastPos[0] + 2e-6]
	ylim = [lastPos[1] - 2e-6, lastPos[1] + 2e-6]

	ax[1].scatter(depositDF['x'][:],depositDF['z'][:],s=10)
	ax[1].scatter(lastPos[0],lastPos[1],label="Particle Center")
	ax[1].legend()
	for axes in ax:
		axes.grid(True)
		axes.set_xlabel('X ($\mu m$)',fontsize=15)
		axes.set_ylabel('Z ($\mu m$)',fontsize=15)
	#ax[1].set_aspect("equal","box")
	ax[1].set_xlim(xlim)
	ax[1].set_ylim(ylim)
	ax[1].set_title("Final orientation of particle")

	#Labels & Legend	
	fig.suptitle(f"Particle simulation",fontsize=20)

#Main() below:

resolution,nDeposits,generateData = parseArgs()
if (generateData == "true"):
	generateNewData()
else:
	print("No new data will be generated, starting plotting... \n")
	
tic = time.time()
print("Starting plotting...")
x,y,z= loadData()




#Poor programming below, please excuse the spaghetti
#Saves the figure in the figures directory, then goes back to generate a new figure and 
#then saves in the figures directory once more.
generateFigure(imageBounds)
os.chdir("..")
os.chdir("figures")

plt.savefig("time_evolution.png")


toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
