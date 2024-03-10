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
	#Compiles and runs the .CPP files
	print(f"Producing {sys.argv[2]} steps \n" )
	subprocess.run(["g++","functions.cpp","particle.cpp","brownian_sim.cpp",
					"-o","sim","-Ofast", "-fopenmp","-funroll-all-loops"])

	subprocess.run(["./sim",nDeposits,sys.argv[3], sys.argv[4],sys.argv[2]])
def compute_MSD(x):
	n = len(x)
	
	MSD = np.zeros((n,1))
	for t in range(n):
		square_sum = 0
		for tau in range(n-t):
			square_sum += (x[tau+t]-x[tau])**2
		MSD[t] = square_sum/(n)
	

	return MSD


def generateFigure():


	fig, ax     = plt.subplots(1, 2, figsize=(10, 7))


	xlim = [- limit,  limit]
	ylim = [- limit*2,  limit]

	#Create circles for particles:
	circle1 = Circle((trajectory[-1,1], trajectory[-1,3]), 2e-6)
	
	circle2 = Circle((trajectory[-1,1], trajectory[-1,2]), 2e-6)
	circle1.set(fill=False, alpha=0.5)
	circle2.set(fill=False, alpha=0.5)
	##############################

	#Generate contour plot of laser profile:

	ax[0].contourf(X,Y,Z,400,alpha=1,zorder=-1)
	ax[1].contourf(X,Y,Z,400,alpha=1,zorder=-1)
	######################################

	#Plot, place, and draw the deposits, circles, trajectory:

	ax[0].plot(trajectory[:,1],trajectory[:,3],label="Trajectory",color='black',zorder=0)
	#ax.plot(particle_2[:,1],particle_2[:,3],label="Trajectory 2",color='blue',zorder=0)
	ax[0].scatter(deposits[-int(nDeposits):,0],deposits[-int(nDeposits):,2],color='red',s=8,alpha=1,zorder=0,label = "Iron oxide")

	ax[0].add_patch(circle1)
	ax[1].add_patch(circle2)

	ax[1].plot(trajectory[:,1],trajectory[:,2],label="Trajectory",color='black',zorder=0)
	ax[1].scatter(deposits[-int(nDeposits):,0],deposits[-int(nDeposits):,1],color='red',s=8,alpha=1,zorder=0,label = "Iron oxide")
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
	plt.savefig("time_evolution.png",facecolor=faceColor)

	os.chdir("..")
	os.chdir("src")

	
	fig, ax     = plt.subplots(3, 1, figsize=(17, 7))

	ax[0].plot(trajectory[:,0],trajectory[:,1],label="X displacement")
	ax[1].plot(trajectory[:,0],trajectory[:,2],label="Y displacement")
	ax[2].plot(trajectory[:,0],trajectory[:,3],label="Z displacement")

	ylabels = ["X","Y","Z"]
	for i in range(3):
		ax[i].set_xlabel("Time (s)")
		ax[i].set_ylabel(ylabels[i])
		ax[i].legend()
		ax[i].grid()

	os.chdir("..")
	os.chdir("figures")
	plt.savefig("displacement.png",facecolor=faceColor)
	plt.clf()
	MSDZ = compute_MSD(trajectory[:,3])
	MSDY = compute_MSD(trajectory[:,2])
	MSDX = compute_MSD(trajectory[:,1])

	n =1
	plt.loglog(trajectory[:-n,0],MSDZ[0:-n]*1e6,label = "Z",linewidth=3)
	plt.loglog(trajectory[:-n,0],MSDX[0:-n]*1e6,label = "X",linewidth=3)
	plt.loglog(trajectory[:-n,0],MSDY[0:-n]*1e6,label = "Y",linewidth=3)
	plt.xlabel("t [s]",fontsize=17)
	plt.ylabel(r"MSD $\langle r^2 ( t) \rangle$",fontsize=17)
	plt.tick_params(axis='both', which='major', labelsize=14)
	plt.tick_params(axis='both', which='minor', labelsize=14)
	plt.grid(True)
	plt.legend(fontsize=15)
	plt.title(f"Mean square displacement for {N-1} timesteps. ",fontsize=20)
	plt.savefig("mean_square_displacement.png",facecolor=faceColor)


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


#generateNewData()

#Trajectory.csv contains time,x,y,z,velocity
trajectory = np.genfromtxt("trajectory.csv",delimiter=',',skip_header=1)
N = len(trajectory)
deposits = np.genfromtxt("deposits.csv",delimiter=',',skip_header=1)
x_positions = trajectory[:,1]
y_positions = trajectory[:,3]
x_deposits = deposits[:, 0]
y_deposits = deposits[:, 2]

print("Computing MSD...")
generateFigure()
print("Done, Creating movie...")

fig,ax= plt.subplots()
sc = ax.scatter([], [], color='b', edgecolor='k', marker='o', facecolor='none')
tic = time.time()
#ani = animation.FuncAnimation(fig, update, frames=N, interval=50,blit=True)
#ani.save('particle_animation.mp4', writer='ffmpeg')

toc= time.time()
print("Done after " + str(round(toc-tic)) + " s")


#plotDeposits()