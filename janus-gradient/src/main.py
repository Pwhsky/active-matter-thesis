import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time 
import subprocess
import sys


#resolution = sys.argv[1] #4000 = 151.648 seconds
#coating    = sys.argv[2]

resolution = "4500"
coating =    "0.5"
iterations =  "1000"

subprocess.run(["g++","main.cpp","-o","sim","-O4"])
subprocess.run(["./sim",resolution,coating,iterations])

#%%
#df = np.genfromtxt('gradient.csv',delimiter=',',skip_header=1)
df = pd.read_csv('gradient.csv')
df.sort_values(by=['y'])


scale = 2e-2
#Extract the data from the DataFrame

x = df['x']
y = df['y']
z = df['z']

#Sum the values to project a 3D image to a 2D image.
Fsum = df.groupby(['x','z'])['gradientValue'].sum().reset_index()


F = Fsum["gradientValue"]

#Normalization
F = (F-F.min())/(F.max()-F.min())



# Create a 3D scatter plot
tic = time.time()

fig = plt.figure()
ax = fig.add_subplot(111)

# Scatter plot with colors based on the field values (F)
# cmap='viridis' is used to map colors to a colormap (you can change it as needed)

sc = ax.scatter(Fsum["x"], Fsum["z"],c=F, cmap='plasma', marker='.',s=130)
#sc = ax.scatter(x,y,z,c=df["gradientValue"], cmap='plasma', marker='.',s=100)

cbar = plt.colorbar(sc)
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_facecolor('black')
ax.set_yticklabels([])
ax.set_xticklabels([])



#ax.set_xlim(-scale,scale)    #For slice representation
#ax.set_ylim(-scale,scale) #For cube representation
#ax.set_ylim(0,0)
#ax.set_zlim(-scale,scale)
# Show the plot
plt. savefig("janus.png")
plt.title("Temperature at time = " + iterations )

toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")

