import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import concurrent.futures
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle
from cython_functions import accelerate_histogram2d, gradient_cython

if len(sys.argv) != 3:
	print("Arguments not given, defaulting to: \n resolution = 200,\n nDeposits = 600 \n")
	resolution = "200" 
	nDeposits = "600"
else:
	resolution 	         = sys.argv[1] 
	nDeposits                = str(sys.argv[2])





print(" Starting simulation...\n")
os.chdir("src")




subprocess.run(["g++","functions.cpp","main.cpp","-o","sim","-Ofast", "-fopenmp" , "-funroll-all-loops"])
subprocess.run(["./sim",resolution,nDeposits])


tic = time.time()
#%%
#df = np.genfromtxt('gradient.csv',delimiter=',',skip_header=1)

df = pd.read_csv("gradient.csv",engine="pyarrow")
df.sort_values(by=['z'])


with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    # Submit each Dask Series for concurrent execution
    futures = [
        executor.submit(dfToNumpy, df['x']),
        executor.submit(dfToNumpy, df['y']),
        executor.submit(dfToNumpy, df['z']),
        executor.submit(dfToNumpy, df['gradientValue'])
    ]

    # Wait for all tasks to complete and retrieve the results
    x = futures[0].result()
    y = futures[1].result()
    z = futures[2].result()
    grad = futures[3].result()


#Sum the slices protruding along the y-axis, kep x and z stationary.

limit = 5e-6

fig = plt.figure()
ax = fig.add_subplot(111)

x_bins = np.linspace(-limit,limit,200)
y_bins = np.linspace(-limit,limit,200)


H, xedges, yedges = histogram2d_cython(x, z, grad, x_bins, y_bins)
#H, xedges, yedges        = np.histogram2d(x, z, bins = [x_bins, y_bins], weights = grad)
#H_counts, xedges, yedges = np.histogram2d(x, z, bins = [x_bins, y_bins]) 
#H = H/H_counts




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
cbar = plt.colorbar()
cbar.set_label(f"$\Delta$T")
plt.imshow(H.T, origin='lower',  cmap='plasma',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])





os.chdir("..")
plt.savefig("gradient.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
subprocess.run(["python3","plotDensity.py"])

