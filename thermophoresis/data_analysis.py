import matplotlib.pyplot as plt
import numpy as np
import os
import time 
import subprocess
import sys
from sklearn.linear_model import LinearRegression

os.chdir("src")
faceColor = (1,0.964706,0.901961)

limit = 15e-6
font = 19
pi = 3.14159
particleRadius = 2e-6


def generateNewData(nDeposits, timesteps, periodicity,bounds):
	print(f"Producing {timesteps} steps \n" )
	#Compiles and runs the .CPP files
	subprocess.run(["g++","functions.cpp","particle.cpp","brownian_sim.cpp",
					"-o","sim","-Ofast", "-fopenmp","-funroll-all-loops","-march=native"])
	subprocess.run(["./sim",str(nDeposits),str(bounds), str(periodicity),str(timesteps)])

	
def compute_MSD(x):
	n = len(x)
	MSD = np.zeros((n,1))
	for t in range(n):
		square_sum = 0
		for tau in range(n-t):
			square_sum += (x[tau+t]-x[tau])**2
		MSD[t] = square_sum/(n-t)
	return MSD

def data_analysis(trajectory):

	#Trajectory.csv contains time,x,y,z,velocity

	#Rescale to be in units of micrometers squared
	MSD = np.asarray([compute_MSD(trajectory[:,1])*1e12, 
				      compute_MSD(trajectory[:,2])*1e12,
				  	  compute_MSD(trajectory[:,3])*1e12	])

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
	x = trajectory[100:-400,0]
	coefficients, covariance = np.polyfit(np.log(x),np.log(np.squeeze(MSD[2][100:-400])),1,cov=True)
	errors                   = np.sqrt(np.diag(covariance))

	#from b = 2*log(v) where v is velocity, we can derive velocity with error bars
	#a and b come from coefficients from numpy polyfit.
	#coefficients[0] = 2
	#print(coefficients)
	#print(errors)

	velocity       = np.exp((coefficients[0])/2)
	velocity_error = np.exp((coefficients[0]+errors[0])/2)-np.exp((coefficients[0]-errors[0])/2)

	#print(f"velocity = {round(velocity,5)} ± {round(velocity_error,5)}" )
	return velocity,velocity_error
	
##############################################################################################


bounds        = 5
timesteps     = 5000
deposit_number = [400,450,475,500,
				 525,550,575,600,
				 625,650,675,700,
				 725,750,775,800,
				 825,850,875,900]
velocities = [2.5136975781808886,  2.445356601463794, 2.458232955184125, 2.6451380562640336, 2.5895883544878666, 2.5376388553024527, 2.608278273190855, 2.643675850756448, 2.688509986448934, 2.6289077882257215, 2.65648510051654, 2.6523037538104965, 2.6434481993678656, 2.629211197357292, 2.709985560952113, 2.6711007696359226, 2.709485467448876, 2.6461686214498408, 2.6399559411670475, 2.6539831410935033]
velocities = [1.4372579662361544,  1.5405952070508544, 1.5704735499707936, 1.756237286631951, 1.868597832306805, 1.3016106800300504, 1.6584314347290847, 1.2199194398675723, 2.1978676522351686, 1.113195449067723, 1.2668201901833174, 2.0522164413705535, 1.8018507971120985, 1.0487169082429464, 1.4106794238208469, 1.9648966070001548, 1.7809671553103275, 1.314269822023115, 1.3913627518311138, 1.3038345247512517]
errors     = [0.003522033740621655, 0.00503520513172484, 0.0047583371822637055, 0.0020076924404301977, 0.001535248075055673, 0.004039533263003214, 0.002492088170126916, 0.0017590080899712035, 0.0025657399853122875, 0.0017927395056758222, 0.003051272116764814, 0.0028058747329615663, 0.001734106070752972, 0.0020798171825058276, 0.0008566913994085823, 0.0016801771113872377, 0.0020286020405806404, 0.00218555576099444, 0.0024901912394086345, 0.0010347779173649485]
errors     = [0.008911047246597636, 0.014675998280480274, 0.0021976301775115115, 0.005650104897457098, 0.009265442420447245, 0.016758057205625176, 0.010648465697403298, 0.006040875405314283, 0.00661019221474568, 0.014477062088285253, 0.009529160392867775, 0.0060108855209364265, 0.008944995425969715, 0.013021525002279155, 0.0073313673925996525, 0.009135094128195353, 0.004582741971437931, 0.016920403945731355, 0.006932246628591265, 0.00564976524541505]

#Perform parameter sweep over laser periodicities
tic = time.time()
#for i in range(len(deposits)):
#	nDeposits = deposits[i]
#	generateNewData(nDeposits, timesteps, 80000*1e-9,bounds)
#	trajectory       = np.genfromtxt("trajectory.csv",delimiter=',',skip_header=3)
#	velocity,error   = data_analysis(trajectory)
	
#	velocities.append(velocity)
#	errors.append(error)
#print(velocities)
#print(errors)

deposits = np.asarray(deposit_number)

r = 30e-9
R = 2e-6
particle_volume = (4*pi*R**3)/3
deposit_volume  = (4*pi*r**3)/3

density_FeO   = 5250  #iron oxide density   in Kg/m³
density_Ag    = 19300 #gold density        in kg/m³
density_SiO   = 2320  #silica oxide density in kg/m³

total_iron   = deposits*deposit_volume*density_FeO
total_gold   = deposits*deposit_volume*density_Ag
total_silica = (particle_volume-(deposits*deposit_volume))*density_SiO

mass_ratios  = np.array(total_iron/(total_silica+total_iron))
velocities = np.array(velocities)
model = LinearRegression()
model.fit(mass_ratios.reshape(-1,1),velocities.flatten())
x_new  = np.linspace(0.0025,0.0075,100)
y_pred = model.predict(x_new.reshape(-1,1))
new_errors = np.std(velocities-model.predict(mass_ratios.reshape(-1,1)))

slope = model.coef_[0]

toc = time.time()
os.chdir("..")
os.chdir("figures")
fig, ax = plt.subplots(figsize=(10, 7))
ax.errorbar(mass_ratios,velocities, errors,capsize=4,fmt='o',label="Average propulsion velocity")
ax.plot(x_new,y_pred,'r-',label=f"Linear fit, k= {slope}")
ax.set_ylabel(r"Velocity [$\mu m /s $]",fontsize=20)
ax.set_xlabel(r"$\frac{m_{\rm Fe}}{m_{\rm tot}}$ ",fontsize=20)
ax.set_title("Propulsion velocity and mass ratio",fontsize=20)
ax.legend()
plt.savefig("vel_mass.png")


print("Finished after " + str(round(toc-tic)) + " s")
##################

