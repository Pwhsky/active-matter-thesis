import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time 
import subprocess
import sys

if len(sys.argv) != 3:
	print("Arguments not given, defaulting to: \n resolution = 100,\n 2D representation \n")
	resolution = "100" #100 is fast, 500 is slow loading in spyder.
	representation = "2";
else:
	resolution 	       = sys.argv[1] #4000 = 151.648 seconds
	representation         = str(sys.argv[2])





print(" Starting simulation...\n")
subprocess.run(["g++","main.cpp","-o","sim","-O4", "-fopenmp"])
subprocess.run(["./sim",resolution,representation])

tic = time.time()
#%%
#df = np.genfromtxt('gradient.csv',delimiter=',',skip_header=1)
df = pd.read_csv('gradient.csv')
df.sort_values(by=['y'])


scale = 5e-6
#Extract the data from the DataFrame

x = df['x']
y = df['y']
z = df['z']


fig = plt.figure()
ax = fig.add_subplot(111)

# Scatter plot with colors based on the field values (F)

#If using full 3D:
#if representation == "3":
#	Fsum = df.groupby(['x','z'])['gradientValue'].sum().reset_index()
#	F = Fsum["gradientValue"]
#	sc = ax.scatter(Fsum["x"],Fsum["z"],c=F, cmap='plasma', marker='.',s=130)
#else:
sc = ax.scatter(df["x"],df["z"],c=df["gradientValue"], cmap='plasma', marker='.',s=130)


ax.set_facecolor('black')
cbar = plt.colorbar(sc,label ="$\Delta T$" )

ax.set_xlabel('X $\mu m$')
ax.set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
ax.set_xticklabels([-20,-10,0,10,20])

ax.set_ylabel('Z $\mu m$')
ax.set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
ax.set_yticklabels([min(z),min(z)/2,0,max(z)/2, max(z)])



plt.title("FeO microparticle temperature gradient")

#ax.set_xlim(-scale,scale)    #For slice representation
#ax.set_ylim(-scale,scale) #For cube representation

# Save the plot
plt.savefig("janus.png")
toc = time.time()
print("Plotting finished after " + str(round(toc-tic)) + " s")

