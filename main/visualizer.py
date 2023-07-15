import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plotMarkerSize=16
plotDPI = 250


nParticles = 6
timeSteps = 1500
#Handler to simulate in C++:
subprocess.run(['g++','brownian-particles.cpp','-o','sim','-O4'])
subprocess.run(['./sim',f' {nParticles}',f' {timeSteps}'])

print("Processing data...")
box_size = 100 * 1e-6

# Read the CSV file into a DataFrame
dfHot = pd.read_csv('hotParticles.csv', header=None, names=['x', 'y'],low_memory=False)
dfCold = pd.read_csv('coldParticles.csv', header=None, names=['x', 'y'],low_memory=False)

# Extract the x and y coordinates from the DataFrame
dfHot['x'] = pd.to_numeric(dfHot['x'],errors='coerce')
dfHot['y'] = pd.to_numeric(dfHot['y'],errors='coerce')


dfCold['x'] = pd.to_numeric(dfCold['x'],errors='coerce')
dfCold['y'] = pd.to_numeric(dfCold['y'],errors='coerce')


# Define the number of timesteps (100 rows per timestep)
#NOTE THIS ASSUMES THAT THERE ARE EQUAL AMMOUNTS OF HOT AND COLD PARTICLES
numTimesteps = dfHot.shape[0]//nParticles
#%%
# Create a figure and axis for the scatter plot
fig, ax = plt.subplots(dpi=plotDPI)

# Create an empty scatter plot
scatterHot = ax.scatter([], [],s=plotMarkerSize,marker='o',color='red',label="Hot Particle")
scatterCold = ax.scatter([], [],s=plotMarkerSize,marker='o',color='blue',label="Cold Particle")
# Set the axis limits
ax.set_xlim([0,box_size]) 
ax.set_ylim([0,box_size])

# Set the axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Function to update the scatter plot for each frame
def update(frame):
    start = frame * nParticles
    end = start + nParticles
    scatterCold.set_offsets(dfCold.iloc[start:end, :].values)
    scatterHot.set_offsets(dfHot.iloc[start:end, :].values)
   
    return scatterHot,scatterCold

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=numTimesteps, interval=1)

# Save the animation as a movie
ani.save('particle_movie.gif', writer='pillow')

# Show the plot (optional)
gif_file = "particle_movie.gif"
subprocess.run(["xdg-open",gif_file])
