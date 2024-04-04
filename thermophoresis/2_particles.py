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
faceColor = (1,0.964706,0.901961)

limit = 15e-6
font = 19
pi = 3.14159
particleRadius = 2e-6

nDeposits          = str(sys.argv[1])

spatialPeriodicity = float(sys.argv[4])*1e-9

def generateNewData():
	print(f"Producing {sys.argv[2]} steps \n" )
	#Compiles and runs the .CPP files
	subprocess.run(["g++","functions.cpp","particle.cpp","brownian_sim.cpp",
					"-o","sim","-Ofast", "-fopenmp","-funroll-all-loops","-march=native"])
	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4],sys.argv[2]])

	
def compute_MSD(x):
	n = len(x)
	MSD = np.zeros((n,1))
	for t in range(n):
		square_sum = 0
		for tau in range(n-t):
			square_sum += (x[tau+t]-x[tau])**2
		MSD[t] = square_sum/(n-t)
	return MSD


def generateFigure():
	makeFigure = True
	if makeFigure:
		fig, ax     = plt.subplots(1, 2, figsize=(10, 7))
		xlim = [- limit,  limit]
		ylim = [- limit*2,  limit]

		#Create circles for particles:
		circle1 = Circle((trajectory_1[-1,1], 
						trajectory_1[-1,3]), 2e-6)
		circle2 = Circle((trajectory_1[-1,1], 
						trajectory_1[-1,2]), 2e-6)


		circle3 = Circle((trajectory_2[-1,1], 
						trajectory_2[-1,3]), 2e-6)
		circle4 = Circle((trajectory_2[-1,1], 
						trajectory_2[-1,2]), 2e-6)


		circle1.set(fill=False, alpha=0.5)
		circle2.set(fill=False, alpha=0.5)
		circle3.set(fill=False, alpha=0.5)
		circle4.set(fill=False, alpha=0.5)
		##############################

		#Generate contour plot of laser profile:

		ax[0].contourf(X,Y,Z,400,alpha=1,zorder=-1)
		ax[1].contourf(X,Y,Z,400,alpha=1,zorder=-1)
		######################################

		#Plot, place, and draw the deposits, circles, trajectory:
		ax[0].add_patch(circle1)
		ax[1].add_patch(circle2)

		ax[0].add_patch(circle3)
		ax[1].add_patch(circle4)

		ax[0].plot(trajectory_1[:,1],
				trajectory_1[:,3],
				label="Trajectory",color='black',zorder=1)
		ax[0].plot(trajectory_2[:,1],
					trajectory_2[:,3],color='black',zorder=1)

		#Plot the last deposits position

		ax[0].scatter(deposits_1[-int(nDeposits):,0],
					  deposits_1[-int(nDeposits):,2],
					color='red',s=8,alpha=1,zorder=0,label = "Iron oxide")

		ax[0].scatter(deposits_2[-int(nDeposits):,0],
					  deposits_2[-int(nDeposits):,2],
					color='red',s=8,alpha=1,zorder=0)



		ax[1].plot(trajectory_1[:,1],
				trajectory_1[:,2],
				label="Trajectory",color='black',zorder=1)
		ax[1].plot(trajectory_2[:,1],
					trajectory_2[:,2],color='black',zorder=1)

		ax[1].scatter(deposits_1[-int(nDeposits):,0],
					deposits_1[-int(nDeposits):,1],
					color='red',s=8,alpha=1,zorder=0,label = "Iron oxide")
		ax[1].scatter(deposits_2[-int(nDeposits):,0],
					deposits_2[-int(nDeposits):,1],
					color='red',s=8,alpha=1,zorder=0)
		#ax[1].scatter(deposits_2[-int(nDeposits):,0],
		#			deposits_2[-int(nDeposits):,1],
		#			color='red',s=8,alpha=1,zorder=0)

		###########################################################
		
		ax[0].set_ylabel('Z (m)',fontsize=16)
		ax[1].set_ylabel('Y (m)',fontsize=16)
		fig.suptitle(f"Particle trajectory for Λ = {round(spatialPeriodicity*1000000)} µm",fontsize=20)
		for i in range(2):
			ax[i].set_xlabel('X (m)',fontsize=16)
			ax[i].axis("equal")
			ax[i].set_xlim(xlim)
			ax[i].set_ylim(ylim)
			ax[i].legend(fontsize=15)
			ax[i].tick_params(axis='both',which='major',labelsize=15)

		#Save figure in figures directory:
		os.chdir("..")
		os.chdir("figures")
		plt.savefig(f"time_evolution_{round(spatialPeriodicity/1000)}_microns.png")#,facecolor=faceColor)
		os.chdir("..")
		os.chdir("src")

##############################################################################################
def update(frame):
	ax.clear()
	ax.set_title('Particle positions at t = {:.2f}'.format(frame))
	ax.set_xlabel('X')
	ax.set_ylabel('Z')
	ax.axis('equal')
	
	crl1 = Circle((trajectory_1[frame,1],trajectory_1[frame,3]),
						radius=2e-6,edgecolor='b',facecolor='none',zorder= 1 )

	crl2 = Circle((trajectory_2[frame,1],trajectory_2[frame,3]),
						radius=2e-6,edgecolor='b',facecolor='none',zorder = 1 )

	ax.add_patch(crl1)
	ax.add_patch(crl2)

	start_index = (frame) * int(nDeposits)
	end_index   = (frame + 1) * int(nDeposits)

	scatter_1      = ax.scatter(x_deposits_1[start_index:end_index], 
							y_deposits_1[start_index:end_index], 
							color='red', marker='o',zorder=0)

	scatter_2      = ax.scatter(x_deposits_2[start_index:end_index], 
							y_deposits_2[start_index:end_index], 
							color='red', marker='o',zorder=0)

	ax.contourf(X,Y,Z,25,alpha=1,zorder=-1)
	ax.set_xlim(-limit, limit)
	ax.set_ylim(-limit*2,limit)



x   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
y   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
X,Y = np.meshgrid(x,y)	
Z = (1+ (np.cos(2*pi*((X))/(spatialPeriodicity))))


################# C++ SIMULATION, YOU CAN REUSE OLD DATA TO SAVE TIME
generateNewData()
##################

#Trajectory.csv contains time,x,y,z,velocity!!!!  
trajectory_1 = np.genfromtxt("trajectory_1.csv",delimiter=',',skip_header=1)
trajectory_2 = np.genfromtxt("trajectory_2.csv",delimiter=',',skip_header=1)

deposits_1   = np.genfromtxt("deposits_1.csv",delimiter=',',skip_header=1)
deposits_2   = np.genfromtxt("deposits_2.csv",delimiter=',',skip_header=1)

#Concatenate deposit data
N =len(trajectory_1)

x_deposits_1  = deposits_1[:, 0]
y_deposits_1  = deposits_1[:, 2]

x_deposits_2  = deposits_2[:, 0]
y_deposits_2  = deposits_2[:, 2]


################ GRAPHS AND FIGURES
generateFigure()
################
print("Creating movie...")

fig,ax = plt.subplots()
sc     = ax.scatter([], [], color='b', edgecolor='k', marker='o', facecolor='none',zorder=0)
tic    = time.time()

################ MOVIE, TAKES A LONG TIME
ani = animation.FuncAnimation(fig, update, frames=N, interval=20,blit=True)
ani.save('particle_animation.mp4', writer='ffmpeg')
################

toc = time.time()

print("Done after " + str(round(toc-tic)) + " s")
