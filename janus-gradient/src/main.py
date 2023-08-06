import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time 
import subprocess
subprocess.run(["g++","main.cpp","-o","sim","-O4"])
subprocess.run(["./sim","4000"])

df = np.genfromtxt('gradient.csv',delimiter=',',skip_header=1)
dfcap = np.genfromtxt('cap.csv',delimiter=',',skip_header=1)

scale = 2e-2
#Extract the data from the DataFrame

x = df[:,0]
y = df[:,1]
z = df[:,2]
F = df[:,3]

F = (F-F.min())/(F.max()-F.min())

#%%
# Create a 3D scatter plot
tic = time.time()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with colors based on the field values (F)
# cmap='viridis' is used to map colors to a colormap (you can change it as needed)

sc = ax.scatter(x, y, z,c=F, cmap='viridis', marker='.',s=130)

cbar = plt.colorbar(sc)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_yticklabels([])
ax.set_xticklabels([])

ax.set_xlim(-scale,scale)    #For slice representation
ax.set_zlim(-scale,scale) #For cube representation
ax.set_ylim(-scale,scale)
#ax.set_zlim(-scale,scale)
ax.grid(True)
# Show the plot
plt.show()
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")