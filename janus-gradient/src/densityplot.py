import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time 
import subprocess
import sys
#TODO: MAKE A HEATMAP DENSITY PLOT THAT LOOKS SICK AF
df = pd.read_csv('deposits.csv')
df.sort_values(by=['y'])


scale = 5e-6
#Extract the data from the DataFrame

x = df['x']
y = df['y']
z = df['z']

fig = plt.figure()
ax = plt.axes(projection="3d")

# Scatter plot with colors based on the field values (F)

#If using full 3D:
#if representation == "3":
#	Fsum = df.groupby(['x','z'])['gradientValue'].sum().reset_index()
#	F = Fsum["gradientValue"]
#	sc = ax.scatter(Fsum["x"],Fsum["z"],c=F, cmap='plasma', marker='.',s=130)

sc = ax.scatter3D(df["x"],df["z"],df["y"], marker='.',s=10)






plt.title("deposit density")

#ax.set_xlim(-scale,scale)    #For slice representation
#ax.set_ylim(-scale,scale) #For cube representation

# Save the plot
plt.savefig("density.jpg")
