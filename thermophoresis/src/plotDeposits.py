import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from matplotlib.patches import Circle




fig, ax     = plt.subplots(1,2 , figsize=(15, 7))
df = pd.read_csv('deposits.csv')
x = df['x']
y = df['y']
z = df['z']
for axes in ax: 
    axes.grid(True)
    axes.tick_params(labelsize=18)
    axes.add_patch(Circle((0, 0), 2e-6,linestyle="--",fill=False))
    #axes.set_xlim([min(df['x'])-1.5e-6,max(df['x'])+1.5e-6  ])
    #axes.set_ylim([min(df['y'])-1.5e-6,max(df['y'])+1.5e-6  ])
    axes.axis('equal')



ax[0].scatter(df['x'],df['z'])
ax[0].set_xlabel("X",fontsize=15)
ax[0].set_ylabel("Z",fontsize=15)
ax[0].set_xlim(-2.5e-6,2.5e-6)
ax[0].set_ylim(-2.5e-6,2.5e-6)

ax[1].scatter(df['y'],df['z'])
ax[1].set_xlabel("Y",fontsize=15)
ax[1].set_ylabel("Z",fontsize=15)
ax[1].set_xlim(-2.5e-6,2.5e-6)
ax[1].set_ylim(-2.5e-6,2.5e-6)
fig.suptitle("Positions of deposits",fontsize=25)

os.chdir("..")
os.chdir("figures")
plt.savefig("deposits.png")

#os.chdir("..")
#subprocess.run(["python3","plotDensity.py"])
