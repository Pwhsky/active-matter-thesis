import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from matplotlib.patches import Circle
import os
#TODO: MAKE A HEATMAP DENSITY PLOT THAT LOOKS SICK AF
os.chdir("src")
df = pd.read_csv('deposits.csv')



xlabel = 'X [m]'
ylabel = 'Y [m]'
zlabel = 'Z [m]'

circle1 = Circle((0, 0), 2e-6)
circle1.set(fill=False, linestyle='--', alpha=0.2)
circle2 = Circle((0, 0), 2e-6)
circle2.set(fill=False, linestyle='--', alpha=0.2)


fig, ax = plt.subplots(1, 2, figsize=(10, 5))
for axis in ax:
	axis.grid(True)
	axis.axis('equal')
	axis.set_ylabel(zlabel)
	

ax[0].set_xlabel(xlabel)
ax[1].set_xlabel(ylabel)


ax[0].scatter(df['x'], df['z'], label='_nolegend_',s=10)
ax[1].scatter(df['y'], df['z'], label='_nolegend_',s=10)

ax[0].add_patch(circle1)
ax[1].add_patch(circle2)



# Add a legend for the circles
fig.legend([circle1, circle2], ["Particle boundary"], loc='lower left')
fig.suptitle(f"Distribution for {len(df)} deposits")
# Save the figure as an image
os.chdir("..")
plt.savefig("density.jpg")
