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

##IMPORTANT:
#The units on the colorbar are kelvin divided by the step-size of the finite differences.
#When resolution = 200 (for most simulations) that equals 25 nanometers!
#Convert accordingly.

font = 17
pi = 3.14159
particleRadius = 2e-6


imageBounds 	   = float(sys.argv[3])*1e-6
spatialPeriodicity = float(sys.argv[4])*1e-9

periodicity = float(sys.argv[4])/1000

def dfToNumpy(column):
    return column.to_numpy()

def parseArgs():
	generateData    		   = str(sys.argv[2])
	resolution 	               = 300
	nDeposits                  = str(sys.argv[1])
	return resolution,nDeposits,generateData
	
def generateNewData():
	subprocess.run(["g++","functions.cpp","particle.cpp","compute_gradient.cpp","-o","sim","-O3", "-fopenmp" , "-funroll-all-loops"])
	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4]])

def loadData():
	df = pd.read_csv("gradient.csv",engine="pyarrow")
	depositDF = pd.read_csv('deposits.csv')
	df.sort_values(by=['z'])

	with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    		#Parallel data loading
    		futures = [
        	executor.submit(dfToNumpy, df['x']),
        	executor.submit(dfToNumpy, df['y']),
      		executor.submit(dfToNumpy, df['z']),
       		executor.submit(dfToNumpy, df['gradX']),
			executor.submit(dfToNumpy, df['gradZ']),
       		]
	x    = futures[0].result()
	y    = futures[1].result()
	z    = futures[2].result()
	gradX = futures[3].result()
	gradZ = futures[4].result()
	return x,y,z,gradX,gradZ,depositDF
	
def generateLaserProfile(spatialPeriodicity): #Generates gaussian laser profile
	x   = np.linspace(-imageBounds*2,imageBounds*2, 200)
	y   = np.linspace(-imageBounds*2,imageBounds*2, 200)
	X,Y = np.meshgrid(x,y)	
	Z = (1+ (np.cos(2*pi*X/spatialPeriodicity)))/2
	return X,Y,Z
	
def generateFigure(imageBounds):
	fig, ax = plt.subplots(1, 3, figsize=(21, 6))
	axisTitles = [f"$∇T$ for {nDeposits} deposits",f"Position of {nDeposits} deposits",f"Laser intensity for $\Lambda$ = {periodicity} μm"] 
	axisLabelsX = ['X ($\mu m$)','X ($\mu m$)','X ($\mu m$)']
	axisLabelsY = ['Z ($\mu m$)','Z ($\mu m$)','Y ($\mu m$)']

	df = pd.read_csv("positions.csv",engine="pyarrow")

	circle1 = Circle((df['x'][0], df['z'][0]), particleRadius)
	circle1.set(fill=False, linestyle='--', alpha=0.2)
	circle2 = Circle((df['x'][1],df['z'][1]), particleRadius)
	circle2.set(fill=False,linestyle='--',alpha=0.2)
	circles	    = [circle1,circle2]
	
	index = 0
	for axis in ax:
		axis.set_title(axisTitles[index],fontsize=font)
		axis.axis('equal')
		axis.set_xlim(-imageBounds,imageBounds)
		axis.set_ylim(-imageBounds,imageBounds)

		axis.set_xlabel(axisLabelsX[index],fontsize=font-3)
		axis.set_ylabel(axisLabelsY[index],fontsize=font-3)
		axis.set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
		axis.set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
		axis.set_xticklabels([round(min(x)*1e6,1),round(min(x)/2*1e6,1),0,round(max(x)/2*1e6,1), round(max(x)*1e6,1)],fontsize=font)
		axis.set_yticklabels([round(min(z)*1e6,1),round(min(z)/2*1e6,1),0,round(max(z)/2*1e6,1), round(max(z)*1e6,1)],fontsize=font)

		if index == 0:
			
	
			#x y z u v w grad are all incrimented in steps of 3, so to subsample(fewer arrows),
			#  your subsampling factor would need to be a multiple of 3, 48 for example.
			subsampling_factor = 114
			#subsampling_factor = 49
			quiver = ax[0].quiver(x[::subsampling_factor],z[::subsampling_factor],
							      w[::subsampling_factor],u[::subsampling_factor],r[::subsampling_factor],
							 
								cmap='plasma',
								pivot = 'tip'
								)
		


			cbar = plt.colorbar(quiver,ax=ax[0])
			cbar.set_label(f"K/μm")
			axis.add_patch(circles[0])
			axis.add_patch(circles[1])
			ax[0].set_facecolor('white')
			
		if index == 1:

			ax[1].grid(True)	
			ax[1].scatter(depositDF['x'], depositDF['z'], label='_nolegend_',s=10)
			
			
			ax[1] = fig.add_subplot(projection='3d')
			ax[1].view_init(azim=90, elev=10)
			ax[1].set_box_aspect([1, 1, 1])
			ax[1].scatter(depositDF['x'], depositDF['y'], depositDF['z'], label='_nolegend_',s=10)
			ax[1].set_xlim(-imageBounds,imageBounds)
			ax[1].set_ylim(-imageBounds,imageBounds)
			ax[1].set_zlim(-imageBounds,imageBounds)
			ax[1].set_xlabel(axisLabelsX[index])
			ax[1].set_ylabel(axisLabelsY[index])
			#ax[1].set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
			#ax[1].set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
			#ax[1].set_zticks([min(x),min(x)/2,0,max(x)/2,max(x)])
			#ax[1].set_xticklabels([round(min(x)*1e6),round(min(x)/2*1e6),0,round(max(x)/2*1e6), round(max(x)*1e6)],fontsize=15)
			#ax[1].set_yticklabels([round(min(z)*1e6),round(min(z)/2*1e6),0,round(max(z)/2*1e6), round(max(z)*1e6)],fontsize=15)
			#ax[1].set_zticklabels([round(min(z)*1e6),round(min(z)/2*1e6),0,round(max(z)/2*1e6), round(max(z)*1e6)],fontsize=15)


		if index == 2:
			X,Y,Z = generateLaserProfile(spatialPeriodicity)
			laserImage = ax[2].contourf(X,Y,Z,50)
			cbar2 = plt.colorbar(laserImage,ax=ax[2])
			cbar2.set_label(" I(x) / $\mathrm{I}_0$")
			
		
		index+=1
	
	#Labels & Legend	
	
	fig.suptitle(f"Silica microparticle temperature gradient",fontsize=20)


#Main() below:
os.chdir("src")
resolution,nDeposits,generateData = parseArgs()

if (generateData == "true"):
	print("Generating new data...\n")
	generateNewData()
else:
	print("No new data will be generated, starting plotting... \n")
	
tic = time.time()
print("Starting plotting...")
x,y,z,gradX,gradZ,depositDF = loadData()

#Scale the sizes with r
r     = np.sqrt(gradX**2 + gradZ**2)
#theta = np.arctan2(gradX,gradZ)

w = gradX
u = gradZ

#unused stuff:

#w     = r*np.sin(theta+pi/2)
#u     = r*np.cos(theta+pi/2)
#tot   = gradX + gradZ

generateFigure(imageBounds)

os.chdir("..")
os.chdir("figures")
plt.savefig("quiver.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
#os.chdir("..")
#subprocess.run(["python3","plotDeos.chdir("..")nsity.py"])


