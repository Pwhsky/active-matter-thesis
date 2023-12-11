
import numpy as np
import matplotlib.pyplot as plt

pi	= 3.14159265358979323846
Dt  = 2e-13
Dr  = 1
v = 3e-6
font = 18
timesteps = 500

positions = np.empty([2,timesteps])
phi       = np.empty([timesteps])
dt = 0.01
positions[0][0] = 0.0
positions[1][0] = 0.0
phi[0] = 0
for t in range(1,timesteps):
    Wx = np.random.normal(loc=0,scale=1)
    Wy = np.random.normal(loc=0,scale=1)
    Wphi = np.random.normal(loc=0,scale=1)
    
    phi[t] = phi[t-1]+ np.sqrt(2*Dr)*Wphi*np.sqrt(dt)
    positions[0][t] = positions[0][t-1] + v*np.cos(phi[t])*dt + np.sqrt(2*Dt)*Wx*np.sqrt(dt)
    positions[1][t] = positions[1][t-1] + v*np.sin(phi[t])*dt + np.sqrt(2*Dt)*Wy*np.sqrt(dt)


fig,ax = plt.subplots(1,figsize=(10,10))
factor = 1e5
x= positions[0][:]*factor
z= positions[1][:]*factor
ax.set_xlim(min(positions[0][:]*factor),max(positions[0][:]*factor))
ax.set_ylim(min(positions[1][:]*factor),max(positions[1][:]*factor))
#ax.set_yticks([min(z),min(z)/2,0,max(z)/2,max(z)])
#ax.set_xticks([min(x),min(x)/2,0,max(x)/2,max(x)])
#ax.set_xticklabels([round(min(x)*factor,1),round(min(x)/2*factor,1),0,round(max(x)/2*factor,1), round(max(x)*factor,1)],fontsize=font)
#ax.set_yticklabels([round(min(z)*factor,1),round(min(z)/2*factor,1),0,round(max(z)/2*factor,1), round(max(z)*factor,1)],fontsize=font)
ax.set_xlabel('X ($\mu m$)',fontsize=17)
ax.set_ylabel('Y ($\mu m$)',fontsize=17)
ax.axis('equal')
ax.scatter(0.0,0.0,label="Start")
ax.scatter(x[-1],z[-1],label="End")
plt.xticks(fontsize=font)
plt.yticks(fontsize=font)
ax.plot(x,z,label="Particle trajectory")
ax.legend(fontsize=15)
ax.grid()
plt.savefig("brownian.pdf")