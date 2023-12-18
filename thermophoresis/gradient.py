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
	depositDF = pd.read_csv('deposits.csv')
	df.sort_values(by=['z'])

	with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    		#Parallel data loading
    		futures = [
        	executor.submit(pandasToNumpy, df['x']),
        	executor.submit(pandasToNumpy, df['y']),
      		executor.submit(pandasToNumpy, df['z']),
       		executor.submit(pandasToNumpy, df['gradX']),
			executor.submit(pandasToNumpy, df['gradZ'])
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
	Z   = (1+ (np.cos(2*pi*X/spatialPeriodicity)))/2
	return X,Y,Z
	
def generateFigureGradient(imageBounds):
	fig, ax     = plt.subplots(1, 3, figsize=(26, 6))
	axisTitles  = [f"∇xT",f"∇zT ",f"Positions for {nDeposits} deposits"] 
	axisLabelsX = ['X ($\mu m$)','X ($\mu m$)','X ($\mu m$)']
	axisLabelsY = ['Z ($\mu m$)','Z ($\mu m$)','Y ($\mu m$)']
	circles	    = [circle1,circle2]
	
	index = 0
	for axis in ax:
		axis.set_title(axisTitles[index],fontsize=font)
		axis.axis('equal')
		axis.set_xlim(-imageBounds,imageBounds)
		axis.set_ylim(-imageBounds,imageBounds)
		
		axis.set_xlabel(axisLabelsX[index],fontsize=font-1)
		axis.set_ylabel(axisLabelsY[index],fontsize=font-1)
		
		axis.set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
		axis.set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
		
		axis.set_xticklabels([round(min(x)*1e6,1),round(min(x)/2*1e6,1),0,round(max(x)/2*1e6,1), round(max(x)*1e6,1)],fontsize=font)
		axis.set_yticklabels([round(min(z)*1e6,1),round(min(z)/2*1e6,1),0,round(max(z)/2*1e6,1), round(max(z)*1e6,1)],fontsize=font)
		if index == 0:
			
			im = ax[0].imshow(H_x.T, origin='lower',  cmap='plasma',
           			 extent=[xedges_x[0], xedges_x[-1], yedges_x[0], yedges_x[-1]])
			cbar = plt.colorbar(im,ax=ax[0])
			cbar.set_label(f"K/μm",fontsize=font-5)
			#axis.add_patch(circles[index])
			
		if index == 1:
			im = ax[1].imshow(H_z.T, origin='lower',  cmap='plasma',
           			 extent=[xedges_z[0], xedges_z[-1], yedges_z[0], yedges_z[-1]])
			cbar = plt.colorbar(im,ax=ax[1])
			cbar.set_label(f"K/μm",fontsize=font-5)

		if index == 2:
			ax[2].grid(True)	
			ax[2].scatter(depositDF['x'], depositDF['z'], label='_nolegend_',s=10)
			axis.add_patch(circles[1])
				
			#X,Y,Z = generateLaserProfile(spatialPeriodicity)
			#laserImage = ax[2].contourf(X,Y,Z,50)
			#cbar2 = plt.colorbar(laserImage,ax=ax[2])
			#cbar2.set_label(" I(x) / $\mathrm{I}_0$")
		index+=1
	
	#Labels & Legend	
	fig.suptitle(f"Silica microparticle temperature gradient",fontsize=20)

def generateFigureQuiver(imageBounds):
	fig, ax = plt.subplots(1, 3, figsize=(26, 6))
	axisTitles = [f"$∇T$ for {nDeposits} deposits",f"Position of {nDeposits} deposits",f"Laser intensity for $\Lambda$ = {periodicity} μm"] 
	axisLabelsX = ['X ($\mu m$)','X ($\mu m$)','X ($\mu m$)']
	axisLabelsY = ['Z ($\mu m$)','Z ($\mu m$)','Y ($\mu m$)']

	df = pd.read_csv("positions.csv",engine="pyarrow")

	circle1 = Circle((df['x'][0], df['z'][0]), particleRadius)
	circle1.set(fill=False, linestyle='--', alpha=0.2)
	#circle2 = Circle((df['x'][1],df['z'][1]), particleRadius)
	#circle2.set(fill=False,linestyle='--',alpha=0.2)
	circles	    = [circle1]
	
	index = 0
	for axis in ax:
		axis.set_title(axisTitles[index],fontsize=font)
		axis.axis('equal')
		axis.set_xlim(-imageBounds,imageBounds)
		axis.set_ylim(-imageBounds,imageBounds)

		axis.set_xlabel(axisLabelsX[index],fontsize=font-1)
		axis.set_ylabel(axisLabelsY[index],fontsize=font-1)
		axis.set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
		axis.set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
		axis.set_xticklabels([round(min(x)*1e6,1),round(min(x)/2*1e6,1),0,round(max(x)/2*1e6,1), round(max(x)*1e6,1)],fontsize=font)
		axis.set_yticklabels([round(min(z)*1e6,1),round(min(z)/2*1e6,1),0,round(max(z)/2*1e6,1), round(max(z)*1e6,1)],fontsize=font)

		if index == 0:
			subsampling_factor = 114
			#subsampling_factor = 49
			quiver = ax[0].quiver(x[::subsampling_factor],z[::subsampling_factor],
							      w[::subsampling_factor],u[::subsampling_factor],r[::subsampling_factor],cmap='plasma',
																										  pivot = 'tip')
		

			cbar = plt.colorbar(quiver,ax=ax[0])
			cbar.set_label(f"K/μm",fontsize=font-5)
			axis.add_patch(circles[0])
			#axis.add_patch(circles[1])
			ax[0].set_facecolor('white')
			
		if index == 1:

			ax[1].grid(True)	
			ax[1].scatter(depositDF['x'], depositDF['z'], label='_nolegend_',s=10)
			
			
			#ax[1] = fig.add_subplot(projection='3d')
			#ax[1].view_init(azim=90, elev=10)
			#ax[1].set_box_aspect([1, 1, 1])
			ax[1].scatter(depositDF['x'], depositDF['z'], label='_nolegend_',s=10,color='C0')
			ax[1].set_xlim(-imageBounds,imageBounds)
			ax[1].set_ylim(-imageBounds,imageBounds)
			#ax[1].set_zlim(-imageBounds,imageBounds)
			ax[1].set_xlabel(axisLabelsX[index])
			ax[1].set_ylabel(axisLabelsY[index])
			circle3 = Circle((df['x'][0], df['z'][0]), particleRadius)
			circle3.set(fill=False, linestyle='--', alpha=0.2)
			ax[1].add_patch(circle3)
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
			cbar2.set_label(" I(x) / $\mathrm{I}_0$",fontsize=font-5)
			
		
		index+=1
	
	#Labels & Legend	
	
	fig.suptitle(f"Silica microparticle temperature gradient",fontsize=20)

#Main() below:

resolution,nDeposits,generateData = parseArgs()
if (generateData == "true"):
	generateNewData()
else:
	print("No new data will be generated, starting plotting... \n")
	
tic = time.time()
print("Starting plotting...")
x,y,z,gradX,gradZ,depositDF = loadData()
x_bins = np.linspace(-imageBounds,imageBounds,200)
y_bins = np.linspace(-imageBounds,imageBounds,200)

H_x, xedges_x, yedges_x = histogram2d_cython(x, z, gradX, x_bins, y_bins)
H_z, xedges_z, yedges_z = histogram2d_cython(x, z, gradZ, x_bins, y_bins)


#Poor programming below, please excuse the spaghetti
#Saves the figure in the figures directory, then goes back to generate a new figure and 
#then saves in the figures directory once more.
os.chdir("..")
os.chdir("figures")
generateFigureGradient(imageBounds)
plt.savefig("temperatureGradient.png")
os.chdir("..")
os.chdir("src")

r = np.sqrt(gradX**2 + gradZ**2)
w = gradX
u = gradZ
generateFigureQuiver(imageBounds)
os.chdir("..")
os.chdir("figures")
plt.savefig("quiver.png")



toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
