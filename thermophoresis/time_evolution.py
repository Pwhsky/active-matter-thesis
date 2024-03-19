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
		circle1 = Circle((trajectory[-1,1], 
						trajectory[-1,3]), 2e-6)
		
		circle2 = Circle((trajectory[-1,1], 
						trajectory[-1,2]), 2e-6)

		circle1.set(fill=False, alpha=0.5)
		circle2.set(fill=False, alpha=0.5)
		##############################

		#Generate contour plot of laser profile:

		ax[0].contourf(X,Y,Z2,400,alpha=1,zorder=-1)
		ax[1].contourf(X,Y,Z1,400,alpha=1,zorder=-1)
		######################################

		#Plot, place, and draw the deposits, circles, trajectory:
		ax[0].add_patch(circle1)
		ax[1].add_patch(circle2)

		ax[0].plot(trajectory[:,1],
				trajectory[:,3],
				label="Trajectory",color='black',zorder=0)
		ax[0].scatter(deposits[-int(nDeposits):,0],
					deposits[-int(nDeposits):,2],
					color='red',s=8,alpha=1,zorder=0,label = "Iron oxide")
		ax[1].plot(trajectory[:,1],
				trajectory[:,2],
				label="Trajectory",color='black',zorder=0)
		ax[1].scatter(deposits[-int(nDeposits):,0],
					deposits[-int(nDeposits):,1],
					color='red',s=8,alpha=1,zorder=0,label = "Iron oxide")
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
		plt.savefig("time_evolution.png")#,facecolor=faceColor)
		os.chdir("..")
		os.chdir("src")
##################################################################################################
	fig, ax     = plt.subplots(3, 1, figsize=(17, 7))
	fig.suptitle("Particle displacement",fontsize=20)
	ylabels = ["X","Y","Z"]

	ax[0].plot(trajectory[:,0],trajectory[:,1],label="X displacement",linewidth=1.5)
	ax[1].plot(trajectory[:,0],trajectory[:,2],label="Y displacement",linewidth=1.5)
	ax[2].plot(trajectory[:,0],trajectory[:,3],label="Z displacement",linewidth=1.5)

	for i in range(3):
		ax[i].set_xlabel("Time (s)",fontsize=18)
		ax[i].set_ylabel(ylabels[i])
		ax[i].legend()
		ax[i].grid()
		ax[i].tick_params(axis='both',which='major',labelsize=15)

	os.chdir("..")
	os.chdir("figures")
	plt.savefig("displacement.png",facecolor=faceColor)
	plt.clf()
######################################################################################################
	fig.suptitle(f"Mean square displacement for Λ = {round(spatialPeriodicity*1000000,1)}µm ",fontsize=24)

	plt.xlabel("t [s]",fontsize=20)
	plt.ylabel(r"MSD(t) [µm²] ",fontsize=20)

	#Rescale to be in units of micrometer
	MSDZ = compute_MSD(trajectory[:,3])*1e12
	MSDY = compute_MSD(trajectory[:,2])*1e12
	MSDX = compute_MSD(trajectory[:,1])*1e12

	plt.loglog(trajectory[2:,0],MSDZ[2:],        label = "Z",linewidth=3)
	plt.loglog(trajectory[2:-100,0],MSDX[2:-100],label = "X",linewidth=3)
	plt.loglog(trajectory[2:-100,0],MSDY[2:-100],label = "Y",linewidth=3)


	
	slope_z = np.asarray( [  np.log(MSDZ[20]/MSDZ[1])  /np.log(20), 
							 np.log(MSDZ[400]/MSDZ[20])/np.log(380/20),
							 np.log(MSDZ[-1]/MSDZ[400])/np.log(len(MSDZ[400:])/380 )  ])

	slope_y = np.asarray( [  np.log(MSDY[20]/MSDY[1])    /np.log(20), 
							 np.log(MSDY[400]/MSDY[20])  /np.log(380/20),
							 np.log(MSDY[-500]/MSDY[400])/np.log((len(MSDY[400:])-500)/380 )  ])

	slope_x = np.asarray( [  np.log(MSDX[20]/MSDX[1])    /np.log(20), 
							 np.log(MSDX[400]/MSDX[20])  /np.log(380/20),
							 np.log(MSDX[-500]/MSDX[400])/np.log((len(MSDX[400:])-500)/380 )  ])

	printSlopes = False
	if printSlopes:
		print("Slopes of Z:")
		print(slope_z)
		print("Slopes of Y:")
		print(slope_y)
		print("Slopes of X:")
		print(slope_x)
	
	num = 1000
	velocities =  [np.divide(MSDX[num:-100],(trajectory[num:-100,0]**2)),
				   np.divide(MSDY[num:-100],(trajectory[num:-100,0]**2)),
				   np.divide(MSDZ[num:-100],(trajectory[num:-100,0]**2))]

	#Unit of v_z is meter per millisecond (dt = 0.01)
	#Rescale to micrometers per second: 
	print("V_x,V_y,V_z")
	for vel in velocities:
		print(np.mean(vel))
	print("err_x, err_y, err_z")	
	for vel in velocities:
		print(np.std(vel)/np.sqrt((len(MSDX[num:-100]))))

	plt.tick_params(axis='both', which='major', labelsize=20)
	plt.tick_params(axis='both', which='minor', labelsize=20)
	plt.grid(True)
	plt.legend(fontsize=15)
	
	plt.savefig("mean_square_displacement.png")
##############################################################################################
def update(frame):
	ax.clear()
	ax.set_title('Particle positions at t = {:.2f}'.format(frame))
	ax.set_xlabel('X')
	ax.set_ylabel('Z')
	ax.axis('equal')
	
	circle = Circle((x_positions[frame],y_positions[frame]),radius=2e-6,edgecolor='b',facecolor='none' )
	ax.add_patch(circle)

	start_index = (frame) * int(nDeposits)
	end_index   = (frame + 1) * int(nDeposits)

	scatter      = ax.scatter(x_deposits[start_index:end_index], 
					    y_deposits[start_index:end_index], 
						color='red', marker='o',zorder=0)

	ax.contourf(X,Y,Z1,25,alpha=1,zorder=-1)
	ax.set_xlim(-limit, limit)
	ax.set_ylim(-limit*2,limit)
	return circle,scatter

x   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
y   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
X,Y = np.meshgrid(x,y)	
Z1 = (1+ (np.cos(2*pi*((X))/(spatialPeriodicity))))
Z2 = (1+ (np.cos(2*pi*(np.sqrt(X**2 + Y**2))/(spatialPeriodicity))))

################# C++ SIMULATION, YOU CAN REUSE OLD DATA TO SAVE TIME
generateNewData()
##################

#Trajectory.csv contains time,x,y,z,velocity
trajectory = np.genfromtxt("trajectory.csv",delimiter=',',skip_header=1)
deposits   = np.genfromtxt("deposits.csv",delimiter=',',skip_header=1)
N = len(trajectory)

x_positions = trajectory[:,1]
y_positions = trajectory[:,2]

x_deposits  = deposits[:, 0]
y_deposits  = deposits[:, 1]

print("Computing MSD...")
################ GRAPHS AND FIGURES
generateFigure()
################
print("Creating movie...")

fig,ax = plt.subplots()
sc     = ax.scatter([], [], color='b', edgecolor='k', marker='o', facecolor='none')
tic    = time.time()

################ MOVIE, TAKES A LONG TIME
ani = animation.FuncAnimation(fig, update, frames=N, interval=50,blit=True)
ani.save('particle_animation.mp4', writer='ffmpeg')
################

toc = time.time()

print("Done after " + str(round(toc-tic)) + " s")
