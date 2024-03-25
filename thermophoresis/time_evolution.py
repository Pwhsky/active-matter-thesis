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

		ax[0].contourf(X,Y,Z,400,alpha=1,zorder=-1)
		ax[1].contourf(X,Y,Z,400,alpha=1,zorder=-1)
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
		plt.savefig(f"time_evolution_{round(spatialPeriodicity/1000)}_microns.png")#,facecolor=faceColor)
		os.chdir("..")
		os.chdir("src")
##################################################################################################
	

	fig, ax     = plt.subplots(3, 1, figsize=(17, 7))
	fig.suptitle("Particle displacement",fontsize=20)
	ylabels = ["X","Y","Z"]

	#ax[0].plot(trajectory[:,0],trajectory[:,1],label="X displacement",linewidth=1.5)
	#ax[1].plot(trajectory[:,0],trajectory[:,2],label="Y displacement",linewidth=1.5)
	#ax[2].plot(trajectory[:,0],trajectory[:,3],label="Z displacement",linewidth=1.5)

	#for i in range(3):
	#	ax[i].set_xlabel("Time (s)",fontsize=18)
	#	ax[i].set_ylabel(ylabels[i])
	#	ax[i].legend()
	#	ax[i].grid()
	#	ax[i].tick_params(axis='both',which='major',labelsize=15)

	os.chdir("..")
	os.chdir("figures")
	plt.savefig("displacement.png",facecolor=faceColor)
	plt.clf()
######################################################################################################

	fig.suptitle(f"Mean square displacement for Λ = {round(spatialPeriodicity*1000000,1)} µm ",fontsize=24)
	legend_labels = ["X","Y","Z"]
	plt.xlabel("t [s]",fontsize=20)
	plt.ylabel(r"MSD(t) [µm²] ",fontsize=20)

	#MSD contains the mean square displacement in different axis[MSD_x, MSD_y, MSD_z]

	#Plot the MSD:
	t = trajectory[:,0]
	plt.xlim([t[1],np.max(t)])
	for i in range(3):
		plt.loglog(trajectory[1:,0],MSD[i][1:],label=legend_labels[i],linewidth=3)


	#Method 1 to obtain slope of MSD, divide MSD with t^2
	#the slopes are in 3 regions of the plot, t[0:20] t[20:400] etc 
	slope_x = np.asarray( [  np.log(MSD[0][20]/MSD[0][1])    /np.log(20), 
							 np.log(MSD[0][400]/MSD[0][20])  /np.log(380/20),
							 np.log(MSD[0][-500]/MSD[0][400])/np.log((len(MSD[0][400:])-500)/380 )  ])
	slope_y = np.asarray( [  np.log(MSD[1][20]/MSD[1][1])    /np.log(20), 
							 np.log(MSD[1][400]/MSD[1][20])  /np.log(380/20),
							 np.log(MSD[1][-500]/MSD[1][400])/np.log((len(MSD[1][400:])-500)/380 )  ])
	slope_z = np.asarray( [  np.log(MSD[2][20]/MSD[2][1])    /np.log(20), 
							 np.log(MSD[2][400]/MSD[2][20])  /np.log(380/20),
							 np.log(MSD[2][-1]/MSD[2][400])  /np.log(len(MSD[2][400:])/380 )  ])

	#Method 2, perform logarithm polyfit (with error estimate):
	x = trajectory[100:-1,0]
	coefficients, covariance = np.polyfit(np.log(x),np.log(np.squeeze(MSD[2][100:-1])),1,cov=True)
	errors                   = np.sqrt(np.diag(covariance))

	#from b = 2*log(v) where v is velocity, we can derive velocity with error bars
	#a and b come from coefficients from numpy polyfit.
	#coefficients[0] = 2
	#print(coefficients)
	#print(errors)

	velocity       = np.exp((coefficients[0])/2)
	velocity_error = np.exp((coefficients[0]+errors[0])/2)-np.exp((coefficients[0]-errors[0])/2)

	print(f"velocity = {round(velocity,5)} ± {round(velocity_error,5)}" )
	plt.loglog(t,np.exp(np.polyval(coefficients,np.log(t))),'--',label = 
						 f"exp({round(coefficients[0],2)}log(t) {round(coefficients[1],2)})")


	#plt.loglog(t,np.log(slope_z[2])*t**slope_z[2],'--',label = "fit")
	printSlopes = False

	plt.tick_params(axis='both', which='major', labelsize=20)
	plt.tick_params(axis='both', which='minor', labelsize=20)
	plt.grid(True)
	plt.legend(fontsize=15)
	
	plt.savefig(f"MSD_{round(spatialPeriodicity/1000)}_microns.png")
##############################################################################################
def update(frame):
	if (frame%2 == 0):
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

		ax.contourf(X,Y,Z,25,alpha=1,zorder=-1)
		ax.set_xlim(-limit, limit)
		ax.set_ylim(-limit*2,limit)

	if frame % 1000 == 0:
		print(f"finished {frame}/{sys.argv[3]} frames.")
	return circle,scatter

x   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
y   = np.linspace(-spatialPeriodicity*20,spatialPeriodicity*20, 200)
X,Y = np.meshgrid(x,y)	
Z = (1+ (np.cos(2*pi*((X))/(spatialPeriodicity))))


################# C++ SIMULATION, YOU CAN REUSE OLD DATA TO SAVE TIME
generateNewData()
##################

#Trajectory.csv contains time,x,y,z,velocity
trajectory = np.genfromtxt("trajectory.csv",delimiter=',',skip_header=3)
deposits   = np.genfromtxt("deposits.csv",delimiter=',',skip_header=1)
N = len(trajectory)
x_positions = trajectory[:,1]
y_positions = trajectory[:,2]
x_deposits  = deposits[:, 0]
y_deposits  = deposits[:, 1]

	#Rescale to be in units of micrometers squared
MSD = np.asarray([compute_MSD(trajectory[:,1])*1e12, 
				  compute_MSD(trajectory[:,2])*1e12,
				  compute_MSD(trajectory[:,3])*1e12	])


print("Computing MSD...")
################ GRAPHS AND FIGURES
generateFigure()
################
print("Creating movie...")

fig,ax = plt.subplots()
sc     = ax.scatter([], [], color='b', edgecolor='k', marker='o', facecolor='none')
tic    = time.time()

################ MOVIE, TAKES A LONG TIME
#ani = animation.FuncAnimation(fig, update, frames=N, interval=50,blit=True)
#ani.save('particle_animation.mp4', writer='ffmpeg')
################

toc = time.time()

print("Done after " + str(round(toc-tic)) + " s")
