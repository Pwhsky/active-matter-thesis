import pandas as pd
import time
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#Avoid interrupting the program as it creates memory leaks, which can only be resolved
#by restarting the IDE.



#%%

plotMarkerSize = 30
plotDPI        = 200


nHot       = 20
nCold      = 30
nParticles = nHot + nCold
timeSteps  = 500
#Handler to simulate in C++:
subprocess.run(['g++','brownian-particles.cpp','-o','sim','-O4'])
subprocess.run(['./sim',f' {nHot}',f' {nCold}',f' {timeSteps}'])

print("Processing data...")
box_size = 100 * 1e-6

# Read the CSV file into a DataFrame
with open('hotParticles.csv', 'r') as file:
    dfHot = pd.read_csv(file, header=None, names=['x', 'y'],low_memory=True)

with open('coldParticles.csv', 'r') as file:
    dfCold = pd.read_csv(file, header=None, names=['x', 'y'], low_memory=True)

# Extract the x and y coordinates from the DataFrame
dfHot['x'] = pd.to_numeric(dfHot['x'],errors='coerce')
dfHot['y'] = pd.to_numeric(dfHot['y'],errors='coerce')

dfCold['x'] = pd.to_numeric(dfCold['x'],errors='coerce')
dfCold['y'] = pd.to_numeric(dfCold['y'],errors='coerce')



#%%
# Create a figure and axis for the scatter plot

fig, ax = plt.subplots(dpi=plotDPI)

# Create an empty scatter plot
scatterHot  = ax.scatter([], [],s=plotMarkerSize,marker='o',color='red',label="Hot Particle")
scatterCold = ax.scatter([], [],s=plotMarkerSize,marker='o',color='blue',label="Cold Particle")
# Set the axis limits
ax.set_xlim([0,box_size]) 
ax.set_ylim([0,box_size])
ax.axis('equal')
# Set the axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
# Draw the 1-micrometer scale bar

# Function to update the scatter plot for each frame
def update(frame):
    #Limit the number of frames in the gif to reduce loading time
    if frame %3 == 0 and frame %3 !=1:
    #if frame %2 ==0:
        startCold = frame * nCold
        endCold   = startCold + nCold +1
        startHot  = frame * nHot
        endHot    = startHot + nHot+1

        scatterCold.set_offsets(dfCold.iloc[startCold:endCold, :].values)
        scatterHot.set_offsets(dfHot.iloc[startHot:endHot, :].values)
    return scatterHot,scatterCold

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=timeSteps, interval=4)


print("Encoding movie...")
tic = time.time()
# Save the animation as a movie
ani.save('particle_movie.gif', writer='pillow')

# Show the plot (optional)
gif_file = "particle_movie.gif"
toc = time.time()

print("Encoding complete, opening movie after "+ str(round(toc-tic)) + " seconds")

subprocess.run(["xdg-open",gif_file])

#Close and clear stuff to prevent leaks
plt.close(fig)
del ani

