import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa

import concurrent.futures
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle
from cython_functions import histogram2d_cython, gradient_cython,transpose_cython

#if len(sys.argv) != 3:
#	print("Invalid arguments, defaulting to: \n resolution = 100,\n deposits = 2 \n")
#	resolution = "100" #100 is fast, 500 is slow loading in spyder.
#	deposits = "2";
#else:
#	resolution 	       = sys.argv[1] #4000 = 151.648 seconds
#	deposits        	= str(sys.argv[2])


def dfToNumpy(column):
    return column.to_numpy()



os.chdir("src")
tic = time.time()

df = pd.read_csv('gradient.csv', engine='pyarrow')

# Sort the data by 'z' column
df = df.sort_values(by='z')

# Compute NumPy arrays if necessary

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


limit = 2.2e-6



fig = plt.figure()
ax = fig.add_subplot(111)

x_bins = np.linspace(-limit,limit,300)
y_bins = np.linspace(-limit,limit,300)



#H, xedges, yedges = np.histogram2d(x, z, bins = [x_bins, y_bins], weights = grad)
#H_counts, xedges, yedges = np.histogram2d(x, z, bins = [x_bins, y_bins]) 
#H = H/H_counts

H, xedges, yedges = histogram2d_cython(x, z, grad, x_bins, y_bins)
H = transpose_cython(H)

dx = xedges[1] - xedges[0]
dy = yedges[1] - yedges[0]
#grad_x, grad_y = np.gradient(H, dx, dy)
grad_x, grad_y = gradient_cython(H, dx, dy)


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
plt.title(f" $∇_z$($\Delta$T)  ")



plt.imshow( abs(grad_x*-1), origin='lower',  cmap='plasma',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
cbar = plt.colorbar()
cbar.set_label(f"[T $L^{-1}$]")
os.chdir("..")
os.chdir("figures")
plt.savefig("gradientContour.pdf",format="pdf")

toc = time.time()
os.chdir("..")
print("Plotting finished after " + str(round(toc-tic)) + " s")
subprocess.run(["python3","plotDensity.py"])
print("Program finished successfully :) \n")

