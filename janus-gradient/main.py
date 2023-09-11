import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle

if len(sys.argv) != 3:
	print("Arguments not given, defaulting to: \n resolution = 100,\n 2D representation \n")
	resolution = "100" #100 is fast, 500 is slow loading in spyder.
	representation = "2";
else:
	resolution 	       = sys.argv[1] #4000 = 151.648 seconds
	representation         = str(sys.argv[2])





print(" Starting simulation...\n")
os.chdir("src")
subprocess.run(["g++","functions.cpp","main.cpp","-o","sim","-O4", "-fopenmp"])
subprocess.run(["./sim",resolution,representation])


tic = time.time()
#%%
#df = np.genfromtxt('gradient.csv',delimiter=',',skip_header=1)
df = pd.read_csv('gradient.csv')
df.sort_values(by=['y'])
x = df['x']
y = df['y']
z = df['z']



fig = plt.figure()
ax = fig.add_subplot(111)

x_bins = np.linspace(-3.5e-6,3.5e-6,400)
y_bins = np.linspace(-3.5e-6,3.5e-6,400)

H, xedges, yedges = np.histogram2d(x, z, bins = [x_bins, y_bins], weights = df['gradientValue'])
H_counts, xedges, yedges = np.histogram2d(x, z, bins = [x_bins, y_bins]) 
H = H/H_counts



ax.set_facecolor('black')
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
plt.title("FeO microparticle temperature gradient")
plt.imshow(H.T, origin='lower',  cmap='plasma',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar()



os.chdir("..")
plt.savefig("gradient.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
subprocess.run(["python3","plotDensity.py"])

