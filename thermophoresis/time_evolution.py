import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time 
import subprocess
import sys
from matplotlib.patches import Circle
import matplotlib.animation as animation

os.chdir("src")


limit = 15e-6
font = 19
pi = 3.14159
particleRadius = 2e-6

nDeposits          = str(sys.argv[1])

spatialPeriodicity = float(sys.argv[4])*1e-9

def generateNewData():
	#Compiles and runs the .CPP files
	print(f"Producing {sys.argv[2]} steps \n" )
	subprocess.run(["g++","functions.cpp","particle.cpp","brownian_sim.cpp",
					"-o","sim","-Ofast", "-fopenmp","-funroll-all-loops"])

	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4],sys.argv[2]])

def generateFigure():

	particle_1     = pd.read_csv("particle_1.csv",engine="pyarrow")
	kinematics     = pd.read_csv("velocity_1.csv",engine="pyarrow")

	#particle_2     = pd.read_csv("particle_2.csv",engine="pyarrow")
	deposits       = pd.read_csv('deposits.csv')

	fig, ax     = plt.subplots(1, 2, figsize=(10, 7))
	#ax.set_title(f"Time evolution after {sys.argv[2]} timesteps.")

	xlim = [- limit,  limit]
	ylim = [- limit*2,  limit]

	#Create circles for particles:
	circle1 = Circle((particle_1.iloc[-1]['x'], particle_1.iloc[-1]['z']), 2e-6)
	
	circle2 = Circle((particle_1.iloc[-1]['x'], particle_1.iloc[-1]['y']), 2e-6)
	circle1.set(fill=False, alpha=0.5)
	circle2.set(fill=False, alpha=0.5)
	##############################

	#Generate contour plot of laser profile:

	ax[0].contourf(X,Y,Z,400,alpha=1,zorder=-1)
	ax[1].contourf(X,Y,Z,400,alpha=1,zorder=-1)
	######################################

	#Plot, place, and draw the deposits, circles, trajectory:
	ax[0].plot(particle_1['x'][:],particle_1["z"][:],label="Trajectory 1",color='black',zorder=0)
	#ax.plot(particle_2['x'][:],particle_2["z"][:],label="Trajectory 2",color='blue',zorder=0)
	ax[0].scatter(deposits['x'][-int(nDeposits):],deposits['z'][-int(nDeposits):],color='red',s=8,alpha=1,zorder=0)

	ax[0].add_patch(circle1)
	ax[1].add_patch(circle2)

	ax[1].plot(particle_1['x'][:],particle_1["y"][:],label="Trajectory 1",color='black',zorder=0)
	ax[1].scatter(deposits['x'][-int(nDeposits):],deposits['y'][-int(nDeposits):],color='red',s=8,alpha=1,zorder=0)
	#ax.add_patch(circle2)

	###########################################################
	
	ax[0].set_ylabel('Z (m)',fontsize=15)
	ax[1].set_ylabel('Y (m)',fontsize=15)

	for i in range(2):
		ax[i].set_xlabel('X (m)',fontsize=15)
		ax[i].axis("equal")
		ax[i].set_xlim(xlim)
		ax[i].set_ylim(ylim)
		ax[i].legend()

	#Save figure in figures directory:
	os.chdir("..")
	os.chdir("figures")
	plt.savefig("time_evolution.png")

	os.chdir("..")
	os.chdir("src")

	
	fig, ax     = plt.subplots(3, 1, figsize=(17, 7))

	ax[0].plot(kinematics['t'],kinematics['x'],label="Displacement")
	ax[1].plot(kinematics['t'],kinematics['y'],label="Displacement")
	ax[2].plot(kinematics['t'],kinematics['z'],label="Displacement")
	#ax[1].plot(kinematics['t'],kinematics['v'],label="Velocity")
	#ax[1].plot(kinematics['t'],kinematics['y'])
	#ax[2].plot(kinematics['t'],kinematics['z'])
	
	ylabels = ["X","Y","Z"]
	for i in range(3):
		ax[i].set_xlabel("Time (s)")
		ax[i].set_ylabel(ylabels[i])
		ax[i].legend()
		ax[i].grid()


	os.chdir("..")
	os.chdir("figures")
	plt.savefig("kinematics.png")

	#################################
def update(frame):
	ax.clear()

	ax.set_xlabel('X')
	ax.set_ylabel('Z')
	ax.set_title('Particle positions at t = {:.2f}'.format(frame))
	ax.axis('equal')
	
	circle = Circle((x_positions[frame],y_positions[frame]),radius=2e-6,edgecolor='b',facecolor='none' )
	ax.add_patch(circle)
	start_index = (frame) * int(nDeposits)
	end_index = (frame + 1) * int(nDeposits)
	scatter = ax.scatter(x_deposits[start_index:end_index], y_deposits[start_index:end_index], color='red', marker='o',zorder=0)
	ax.contourf(X,Y,Z,25,alpha=1,zorder=-1)
	
	ax.set_xlim(-limit, limit)
	ax.set_ylim(-limit*2,limit)

	return circle,scatter

x   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
y   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
X,Y = np.meshgrid(x,y)	
Z = (1+ (np.cos(2*pi*X/(spatialPeriodicity))))


generateNewData()

data = np.genfromtxt("particle_1.csv",delimiter=',',skip_header=1)
deposits = np.genfromtxt("deposits.csv",delimiter=',',skip_header=1)

x_positions = data[:, 0]
y_positions = data[:, 2]

x_deposits = deposits[:, 0]
y_deposits = deposits[:, 2]


fig,ax= plt.subplots()
sc = ax.scatter([], [], color='b', edgecolor='k', marker='o', facecolor='none')


print("Creating movie...")
tic = time.time()

ani = animation.FuncAnimation(fig, update, frames=len(data), interval=50,blit=True)
ani.save('particle_animation.mp4', writer='ffmpeg')


toc= time.time()
print("Movie created after " + str(round(toc-tic)) + " s")
generateFigure()

#plotDeposits()