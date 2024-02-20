import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from matplotlib.patches import Circle


fig, ax     = plt.subplots(1, 3, figsize=(13, 7))
df = pd.read_csv('deposits.csv')
x = df['x']
y=df['y']
z=df['z']
ax = fig.add_subplot(projection='3d')
ax.view_init(azim=60, elev=10)
ax.set_box_aspect([1, 1, 1])

ax.set_xlim([min(df['x'])-2e-6,max(df['x'])+1e-6  ])
ax.set_ylim([min(df['y'])-2e-6,max(df['y'])+1e-6  ])
ax.set_zlim([min(df['z'])-2e-6,max(df['z'])+1e-6  ])

ax.scatter3D(df['x'],df['y'],df['z'])

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
os.chdir("..")
os.chdir("figures")
plt.savefig("deposits.png")

#os.chdir("..")
#subprocess.run(["python3","plotDensity.py"])
