import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time 
import subprocess
import sys
from matplotlib.patches import Circle
#TODO: MAKE A HEATMAP DENSITY PLOT THAT LOOKS SICK AF
df = pd.read_csv('deposits.csv')




num_bins = 50

max_distance = 5e-6
bin_width = max_distance/num_bins

rdf_values = np.zeros(num_bins)

distances = np.linalg.norm(df[['x', 'y', 'z']].values[:, np.newaxis, :] - df[['x', 'y', 'z']].values, axis=2)
np.fill_diagonal(distances, np.inf)


xlabel = 'X [m]'
ylabel = 'Y [m]'
zlabel = 'Z [m]'


fig1 = plt.figure(figsize=(8, 6))
ax = fig1.add_subplot(111)
ax.scatter(df['x'], df['z'],label='_nolegend_')
ax.axis('equal')

plt.xlabel(xlabel)
plt.ylabel(zlabel)
plt.title('FeO deposit distribution')
plt.grid(True)
#overlay with particle radius
circle = Circle((0,0), 2e-6)
circle.set(fill=False,linestyle='--',alpha=0.2)
ax.add_patch(circle)

#Legend
plt.legend([ "Particle boundary"],loc='lower left')
plt.savefig("densityxz.jpg")



fig2 = plt.figure(figsize=(8, 6))
ax = fig2.add_subplot(111)
ax.axis('equal')

ax.scatter(df['y'], df['z'],label='_nolegend_')
plt.xlabel(ylabel)
plt.ylabel(zlabel)
plt.title('FeO deposit distribution')
plt.grid(True)
#overlay with particle radius
circle = Circle((0,0), 2e-6)
circle.set(fill=False,linestyle='--',alpha=0.2)
ax.add_patch(circle)
#Legend
plt.legend([ "Particle boundary"],loc='lower left')
plt.savefig("densityyz.jpg")
