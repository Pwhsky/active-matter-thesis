import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle

if len(sys.argv) != 3:
	print("Arguments not given, defaulting to: \n resolution = 200,\n nDeposits = 600 \n")
	resolution = "200" 
	nDeposits = "600"
else:
	resolution 	       = sys.argv[1] 
	nDeposits                = str(sys.argv[2])





print(" Starting simulation...\n")
os.chdir("src")
subprocess.run(["g++","functions.cpp","main.cpp","-o","sim","-O4", "-fopenmp"])
subprocess.run(["./sim",resolution,nDeposits])


tic = time.time()
#%%
#df = np.genfromtxt('gradient.csv',delimiter=',',skip_header=1)
df = pd.read_csv('gradient.csv')
df.sort_values(by=['z'])
x = df['x'].to_numpy()
y = df['y'].to_numpy()
z = df['z'].to_numpy()
grad = df['gradientValue'].to_numpy()


#Sum the slices protruding along the y-axis, kep x and z stationary.


limit = 5e-6

fig = plt.figure()
ax = fig.add_subplot(111)

x_bins = np.linspace(-limit,limit,400)
y_bins = np.linspace(-limit,limit,400)

H, xedges, yedges        = np.histogram2d(x, z, bins = [x_bins, y_bins], weights = grad)
H_counts, xedges, yedges = np.histogram2d(x, z, bins = [x_bins, y_bins]) 
H = H/H_counts




#cbar = plt.colorbar(sc,label ="$\Delta T$" )
ax.set_xlabel('X [m]')
ax.set_ylabel('Z [m]')


#overlay with particle radius
circle = Circle((0,0), 2e-6)
circle.set(fill=False,linestyle='--',alpha=0.2)
ax.add_patch(circle)

#Legend
plt.legend([ "Particle boundary"],loc='lower left')
# Save the plot
plt.title(f"$\Delta$T for {nDeposits} deposits")
ax.set_facecolor('black')
plt.imshow(H.T, origin='lower',  cmap='plasma',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar()




os.chdir("..")
plt.savefig("gradient.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
subprocess.run(["python3","plotDensity.py"])

