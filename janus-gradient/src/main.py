import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time 
import subprocess
import sys
from matplotlib_scalebar.scalebar import ScaleBar

#resolution = sys.argv[1] #4000 = 151.648 seconds
#coating    = sys.argv[2]
resolution = "100" #100 is fast, 500 is slow loading.


subprocess.run(["g++","main.cpp","-o","sim","-O4", "-fopenmp"])
subprocess.run(["./sim",resolution])

#%%
#df = np.genfromtxt('gradient.csv',delimiter=',',skip_header=1)
df = pd.read_csv('gradient.csv')
df.sort_values(by=['y'])


scale = 5e-6
#Extract the data from the DataFrame

x = df['x']
y = df['y']
z = df['z']

#Sum the values to project a 3D image to a 2D image.
Fsum = df.groupby(['x','z'])['gradientValue'].sum().reset_index()


F = Fsum["gradientValue"]

#Normalization
#F = (F-F.min())/(F.max()-F.min())



# Create a 3D scatter plot
tic = time.time()

fig = plt.figure()
ax = fig.add_subplot(111)

# Scatter plot with colors based on the field values (F)
# cmap='viridis' is used to map colors to a colormap (you can change it as needed)


sc = ax.scatter(df["x"],df["z"],c=df["gradientValue"], cmap='plasma', marker='.',s=130)
#sc = ax.scatter(Fsum["x"],Fsum["z"],c=F, cmap='plasma', marker='.',s=130)

#sc = ax.scatter(x,y,z,c=df["gradientValue"], cmap='plasma', marker='.',s=100)
ax.set_facecolor('black')

cbar = plt.colorbar(sc,label ="$\Delta T$" )
ax.set_xlabel('X $\mu m$')

ax.set_ylabel('Z $\mu m$')

ax.set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
ax.set_yticklabels([-20,-10,0,10,20])

ax.set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
ax.set_xticklabels([-20,-10,0,10,20])


plt.title("FeO microparticle temperature gradient")



ax.set_xlim(-scale,scale)    #For slice representation
ax.set_ylim(-scale,scale) #For cube representation

#ax.set_ylim(0,0)
#ax.set_zlim(-scale,scale)
# Show the plot
plt. savefig("janus.png")


toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")

