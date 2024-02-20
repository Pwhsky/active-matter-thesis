
import numpy as np
import matplotlib.pyplot as plt

pi	= 3.14159265358979323846
Dt  = 2e-13
Dr  = 2
v = 3e-6
font = 18
timesteps = 100000

positions = np.empty([2,timesteps])
phi       = np.empty([timesteps])
dt = 0.01

positions[0][0] = 0.0
positions[1][0] = 0.0
phi[0] = 0
for t in range(1,timesteps):
    Wx   = np.random.normal(loc=0,scale=1)
    Wy   = np.random.normal(loc=0,scale=1)
    Wphi = np.random.normal(loc=0,scale=1)
    
    phi[t] = phi[t-1]+ np.sqrt(2*Dr)*Wphi*np.sqrt(dt)
    positions[0][t] = positions[0][t-1] + v*np.cos(phi[t])*dt + np.sqrt(2*Dt)*Wx*np.sqrt(dt)
    positions[1][t] = positions[1][t-1] + v*np.sin(phi[t])*dt + np.sqrt(2*Dt)*Wy*np.sqrt(dt)

#amogus
fig,ax = plt.subplots(1,figsize=(11,11))
factor = 1e5
x= positions[0][:]*factor
z= positions[1][:]*factor
ax.set_xlim(min(positions[0][:]*factor),max(positions[0][:]*factor))
ax.set_ylim(min(positions[1][:]*factor),max(positions[1][:]*factor))
ax.set_xlabel('X ($\mu m$)',fontsize=19)
ax.set_ylabel('Y ($\mu m$)',fontsize=19)
ax.axis('equal')
ax.scatter(0.0,0.0,label="Start",c='r')
ax.scatter(x[-1],z[-1],label="End",c='k')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
ax.plot(x,z,label="Particle trajectory",linewidth=1.5)
ax.legend(fontsize=17)
ax.grid()
plt.savefig("brownian.png")