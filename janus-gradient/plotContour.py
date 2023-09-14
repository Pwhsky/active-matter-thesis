import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle

#if len(sys.argv) != 3:
#	print("Invalid arguments, defaulting to: \n resolution = 100,\n deposits = 2 \n")
#	resolution = "100" #100 is fast, 500 is slow loading in spyder.
#	deposits = "2";
#else:
#	resolution 	       = sys.argv[1] #4000 = 151.648 seconds
#	deposits        	= str(sys.argv[2])





#print(" Starting simulation...\n")
os.chdir("src")
#subprocess.run(["g++","functions.cpp","main.cpp","-o","sim","-O4", "-fopenmp"])
#subprocess.run(["./sim",resolution,deposits])


tic = time.time()
#%%
#df = np.genfromtxt('gradient.csv',delimiter=',',skip_header=1)
df = pd.read_csv('gradient.csv')
df.sort_values(by=['z'])
x = df['x'].to_numpy()
y = df['y'].to_numpy()
z = df['z'].to_numpy()
grad = df['gradientValue'].to_numpy()

limit = 5e-6



fig = plt.figure()
ax = fig.add_subplot(111)

x_bins = np.linspace(-limit,limit,200)
y_bins = np.linspace(-limit,limit,200)



H, xedges, yedges = np.histogram2d(x, z, bins = [x_bins, y_bins], weights = grad)
H_counts, xedges, yedges = np.histogram2d(x, z, bins = [x_bins, y_bins]) 
H = H/H_counts
H = H.T

dx = xedges[1] - xedges[0]
dy = yedges[1] - yedges[0]
grad_x, grad_y = np.gradient(H, dx, dy,edge_order=1)

# Create a meshgrid for the gradient vectors
X_grad, Y_grad = np.meshgrid(xedges[0:-1], yedges[0:-1])




#cbar = plt.colorbar(sc,label ="$\Delta T$" )
ax.set_xlabel('X [m]')
ax.set_ylabel('Z [m]')


#overlay with particle radius
circle = Circle((0,0), 2e-6)
circle.set(fill=False,linestyle='--',alpha=0.2)
ax.add_patch(circle)
ax.axis('equal')
ax.set_facecolor('black')

#Legend
plt.legend([ "Particle boundary"],loc='lower left')
# Save the plot
plt.title(f"gradient of temperature gradient $\Delta$T ")
os.chdir("..")

plt.imshow(grad_x, origin='lower',  cmap='plasma',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar()

plt.savefig("gradientContour.png")

toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
subprocess.run(["python3","plotDensity.py"])
print("Program finished successfully :) \n")

