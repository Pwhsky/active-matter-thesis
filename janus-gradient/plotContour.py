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

pi = 3.14159
circle1 = Circle((0, 0), 2e-6)
circle1.set(fill=False, linestyle='--', alpha=0.2)
circle2 = Circle((0,0), 2e-6)
circle2.set(fill=False,linestyle='--',alpha=0.2)

imageBounds 	   = 2.3*1e-6
spatialPeriodicity   = 80000*1e-9

def dfToNumpy(column):
    return column.to_numpy()

	

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
       		executor.submit(dfToNumpy, df['gradientValue'])
       		]
	x    = futures[0].result()
	y    = futures[1].result()
	z    = futures[2].result()
	grad = futures[3].result()
	return x,y,z,grad,depositDF
	

def generateFigure(imageBounds):
	fig, axis    = plt.subplots( )
	axisTitles   = [f"Contour plot"] 
	axisLabelsX  = ['X ($\mu m$)']
	axisLabelsY  = ['Z ($\mu m$)']
	circles	    = [circle1]
	index = 0
	axis.set_title(axisTitles[index])
	axis.axis('equal')
	axis.set_xlim(-imageBounds,imageBounds)
	axis.set_ylim(-imageBounds,imageBounds)
		
	axis.set_xlabel(axisLabelsX[index])
	axis.set_ylabel(axisLabelsY[index])
		
	axis.set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
	axis.set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
	
	axis.set_xticklabels([round(min(x)*1e6),round(min(x)/2*1e6),0,round(max(x)/2*1e6), round(max(x)*1e6)])
	axis.set_yticklabels([round(min(z)*1e6),round(min(z)/2*1e6),0,round(max(z)/2*1e6), round(max(z)*1e6)])
	axis.add_patch(circles[index])
	im = axis.imshow(H.T, origin='lower',  cmap='plasma',
           			 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	cbar = plt.colorbar(im,axis)
	cbar.set_label(f"$\Delta$T [K]")
	
	axis.add_patch(circles[index])

	#Labels & Legend	
	fig.legend([ "Particle boundary"],loc='lower left')
	fig.suptitle(f"Silica microparticle temperature gradient")
	
#Main() below:
os.chdir("src")
print("Plotting contours...")

	
tic = time.time()

x,y,z,grad,depositDF = loadData()


x_bins = np.linspace(-imageBounds,imageBounds,200)
y_bins = np.linspace(-imageBounds,imageBounds,200)
H, xedges, yedges = histogram2d_cython(x, z, grad, x_bins, y_bins)
generateFigure(imageBounds)

os.chdir("..")
os.chdir("figures")
plt.savefig("contour.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
#os.chdir("..")
#subprocess.run(["python3","plotDensity.py"])


