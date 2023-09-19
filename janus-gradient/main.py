import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import concurrent.futures
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle
from cython_functions import histogram2d_cython, gradient_cython
imageBounds = 2.2e-6


def dfToNumpy(column):
    return column.to_numpy()

def parseArgs():
	if len(sys.argv) != 4:
		userInput    = input("Would you like to generate new data? true/false \n")
		resolution   = "150" 
		nDeposits    = "600"
		generateData = userInput
	else:
		resolution 	         = sys.argv[1] 
		nDeposits                = str(sys.argv[2])
		generateData		 = sys.argv[3]
	return resolution,nDeposits,generateData
def generateNewData():
	
	print("Generating new data...\n")
	subprocess.run(["g++","functions.cpp","main.cpp","-o","sim","-Ofast", "-fopenmp" , "-funroll-all-loops"])
	subprocess.run(["./sim",resolution,nDeposits])

def loadData():
	df = pd.read_csv("gradient.csv",engine="pyarrow")
	df.sort_values(by=['z'])

	with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    		#Parallel data loading
    		futures = [
        	executor.submit(dfToNumpy, df['x']),
        	executor.submit(dfToNumpy, df['y']),
      		executor.submit(dfToNumpy, df['z']),
       		executor.submit(dfToNumpy, df['gradientValue'])
    	]
	x    = futures[0].result()
	y    = futures[1].result()
	z    = futures[2].result()
	grad = futures[3].result()
	return x,y,z,grad
	
def generateFigure():
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	x_bins = np.linspace(-imageBounds,imageBounds,200)
	y_bins = np.linspace(-imageBounds,imageBounds,200)
	H, xedges, yedges = histogram2d_cython(x, z, grad, x_bins, y_bins)


	#Add circle
	circle = Circle((0,0), 2e-6)
	circle.set(fill=False,linestyle='--',alpha=0.2)
	ax.add_patch(circle)
	
	#Labels & Legend
	ax.set_xlabel('X [m]')
	ax.set_ylabel('Z [m]')	
	plt.legend([ "Particle boundary"],loc='lower left')
	plt.title(f"$\Delta$T for {nDeposits} deposits")
	ax.set_facecolor('black')
	
	#Colorbar & imshow
	plt.imshow(H.T, origin='lower',  cmap='plasma',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	cbar = plt.colorbar()
	cbar.set_label(f"$\Delta$T")

#Main() below:
os.chdir("src")
resolution,nDeposits,generateData = parseArgs()

if (generateData == "true"):
	generateNewData()
else:
	print("No new data will be generated, starting plotting... \n")
	
tic = time.time()

x,y,z,grad = loadData()
generateFigure()

os.chdir("..")
os.chdir("figures")
plt.savefig("gradient.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
os.chdir("..")
subprocess.run(["python3","plotDensity.py"])


