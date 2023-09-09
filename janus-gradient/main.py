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


# Scatter plot with colors based on the field values (F)

#If using full 3D:
#if representation == "3":
#	Fsum = df.groupby(['x','z'])['gradientValue'].sum().reset_index()
#	F = Fsum["gradientValue"]
#	sc = ax.scatter(Fsum["x"],Fsum["z"],c=F, cmap='plasma', marker='.',s=130)
fig = plt.figure()
ax = fig.add_subplot(111)


sc = ax.scatter(df["x"],df["z"],c=df["gradientValue"], cmap='plasma', marker='.',s=130,label='_nolegend_')


ax.set_facecolor('black')
cbar = plt.colorbar(sc,label ="$\Delta T$" )

ax.set_xlabel('X $\mu m$')
ax.set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
ax.set_xticklabels([round(min(x)*1e6),round(min(x)/2*1e6),0,round(max(x)/2*1e6), round(max(x)*1e6)])

ax.set_ylabel('Z $\mu m$')
ax.set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
ax.set_yticklabels([round(min(z)*1e6),round(min(z)/2*1e6),0,round(max(z)/2*1e6), round(max(z)*1e6)])

#overlay with particle radius
circle = Circle((0,0), 2e-6)
circle.set(fill=False,linestyle='--',alpha=0.2)
ax.add_patch(circle)


#Legend
plt.legend([ "Particle boundary"],loc='lower left')

# Save the plot
plt.title("FeO microparticle temperature gradient")
os.chdir("..")
plt.savefig("gradient.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")
subprocess.run(["python3","plotDensity.py"])
