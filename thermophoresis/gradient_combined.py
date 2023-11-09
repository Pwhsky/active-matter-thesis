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


pi = 3.14159
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
	
def set_equal_axes_limits(ax, image_bounds):
    ax.axis('equal')
    ax.set_xlim(-image_bounds, image_bounds)
    ax.set_ylim(-image_bounds, image_bounds)

def set_labels_and_ticks(ax, x, z):
    ax.set_xlabel('X ($\mu m$)',fontsize=16)
    ax.set_ylabel('Z ($\mu m$)',fontsize=16)
    ax.set_yticks([min(z), min(z) / 2, 0, max(z) / 2, max(z)])
    ax.set_xticks([min(x), min(x) / 2, 0, max(x) / 2, max(x)])
    ax.set_xticklabels([round(val * 1e6, 1) for val in [min(x), min(x) / 2, 0, max(x) / 2, max(x)]], fontsize=15)
    ax.set_yticklabels([round(val * 1e6, 1) for val in [min(z), min(z) / 2, 0, max(z) / 2, max(z)]], fontsize=15)


	
def generateGradient():
	print("Generating new data...\n")
	subprocess.run(["g++","functions.cpp","compute_gradient.cpp","-o","sim","-Ofast", "-fopenmp" , "-funroll-all-loops"])
	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4]])

def generateTemperature():
	print("Generating new data...\n")
	subprocess.run(["g++","functions.cpp","compute_temperature.cpp","-o","sim","-Ofast", "-fopenmp" , "-funroll-all-loops"])
	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4]])

def loadData():
	df1 = pd.read_csv("gradient.csv",engine="pyarrow")
	df2 = pd.read_csv("temperature.csv",engine="pyarrow")
	depositDF = pd.read_csv('deposits.csv')
	#df.sort_values(by=['z'])

	with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    		#Parallel data loading
    		futures = [
        	executor.submit(pandasToNumpy, df1['x']),
        	executor.submit(pandasToNumpy, df1['y']),
      		executor.submit(pandasToNumpy, df1['z']),
       		executor.submit(pandasToNumpy, df1['gradX']),
			executor.submit(pandasToNumpy, df1['gradZ']),
			executor.submit(pandasToNumpy, df2['temperature'])

       		]
	x    = futures[0].result()
	y    = futures[1].result()
	z    = futures[2].result()
	gradX = futures[3].result()
	gradZ = futures[4].result()
	temperature = futures[5].result()
	return x,y,z,gradX,gradZ,temperature,depositDF
	
def generateLaserProfile(spatialPeriodicity): #Generates gaussian laser profile
	x   = np.linspace(-imageBounds*2,imageBounds*2, 200)
	y   = np.linspace(-imageBounds*2,imageBounds*2, 200)
	X,Y = np.meshgrid(x,y)	
	Z   = (1+ (np.cos(2*pi*X/spatialPeriodicity)))/2
	return X,Y,Z
	
def generateFigure(imageBounds):
	fig, axis     = plt.subplots(2, 3, figsize=(30, 15))
	axisTitlesTop  = [f"Temperature increase ΔT",f"∇zT ",f"Illumination profile with spatial periodicity = {periodicity} μm"] 
	axisTitlesBottom  = [f"∇xT",f"∇zT ",f"Tangential flow"] 
	circles	    = [circle1,circle2]

	#Generate tick labels and limits
	for i in range(2):
		for j in range(3):
			if i == 0:
				axis[i,j].set_title(axisTitlesTop[j])
			if i == 1:

				axis[i,j].set_title(axisTitlesBottom[j])
			set_equal_axes_limits(axis[i, j], imageBounds)
			set_labels_and_ticks(axis[i, j], x, z)



	im1 = axis[1,0].imshow(H_x.T, origin='lower',  cmap='plasma',
           	 extent=[xedges_x[0], xedges_x[-1], yedges_x[0], yedges_x[-1]])
	cbar1 = plt.colorbar(im1,ax=axis[0,0])

	cbar1.set_label(f"$\Delta$T [K]")

	
	im2 = axis[1,1].imshow(H_z.T, origin='lower',  cmap='plasma',
    	extent=[xedges_z[0], xedges_z[-1], yedges_z[0], yedges_z[-1]])
	cbar2 = plt.colorbar(im2,ax=axis[1,1])
	cbar2.set_label(f"K/μm")

	
	subsampling_factor = 114
		#subsampling_factor = 49
	im3 = axis[1,2].quiver(x[::subsampling_factor],z[::subsampling_factor],
							  w[::subsampling_factor],u[::subsampling_factor],r[::subsampling_factor],
							 
							cmap='plasma',
							pivot = 'tip'
							)
	cbar3 = plt.colorbar(im3,ax=axis[1,2])
	cbar3.set_label(f"K/μm")
	#axis.add_patch(circles[0])
	#axis.add_patch(circles[1])
	axis[0,2].set_facecolor('white')

	
	im0 = axis[0,0].imshow(H_temp.T,origin='lower',cmap='plasma',
		extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])	
	cbar0 = plt.colorbar(im3,ax=axis[1,0])
	cbar0.set_label(f"K/μm")


	axis[0,0].scatter(depositDF['x'], depositDF['z'], label='Iron Oxide',s=10)
	axis[0,0].legend()
		
	X,Y,Z = generateLaserProfile(spatialPeriodicity)
	laserImage = axis[0,2].contourf(X,Y,Z,50)
	cbar2 = plt.colorbar(laserImage,ax=axis[0,2])
	cbar2.set_label(" I(x) / $\mathrm{I}_0$",fontsize=15)
			
		
	
	
	#Labels & Legend	
	fig.legend([ "Particle boundary"],loc='lower left')
	fig.suptitle(f"Silica microparticle temperature gradient",fontsize=20)


#Main() below:
#Navigate to folder containing c++ program

resolution,nDeposits,generateData = parseArgs()
if (generateData == "true"):
	generateGradient()
	generateTemperature()
else:
	print("No new data will be generated, starting plotting... \n")
	
tic = time.time()
print("Starting plotting...")
x,y,z,gradX,gradZ,temperature,depositDF = loadData()

x_bins = np.linspace(-imageBounds,imageBounds,200)
y_bins = np.linspace(-imageBounds,imageBounds,200)
#Use cython to accelerate histogram generation

H_x, xedges_x, yedges_x = histogram2d_cython(x, z, gradX, x_bins, y_bins)
H_z, xedges_z, yedges_z = histogram2d_cython(x, z, gradZ, x_bins, y_bins)
H_temp, xedges, yedges = histogram2d_cython(x, z, temperature, x_bins, y_bins)
r     = np.sqrt(gradX**2 + gradZ**2)
w = gradX
u = gradZ

generateFigure(imageBounds)

#Save figure to figures directory
os.chdir("..")
os.chdir("figures")
plt.savefig("combined.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
