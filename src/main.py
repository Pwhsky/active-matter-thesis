import pandas as pd
import time
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing
import os

#Avoid interrupting the program as it creates memory leaks, which can only be resolved
#by restarting the IDE.

#The simulation timestep is edited in the CPP source file.

#%%

#Select parameters
plotMarkerSize = 30
plotDPI        = 200

clusterType = "rotator"
timeSteps  = 500

if clusterType == "spinner":
    nHot        = 2
    nCold       = 4
elif clusterType == "rotator":
    nHot        = 2
    nCold       = 3   
elif clusterType == "stator":
    nHot        = 2
    nCold       = 2
elif clusterType == "migrator":
    nHot        = 1
    nCold       = 1
  
nParticles = nHot + nCold

#Handler to start simulation
subprocess.run(['g++','main.cpp','particle.cpp','-o','sim','-O4'])
subprocess.run(['./sim',clusterType,f' {timeSteps}'])
print("Processing data...")




# Read the CSV files into a pandas dataframe
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

#Prepare the scatter plot:
box_size = 100 * 1e-6
fig, ax = plt.subplots(dpi=plotDPI)

scatterHot  = ax.scatter([], [],s=plotMarkerSize,marker='o',color='red',label="Hot Particle")
scatterCold = ax.scatter([], [],s=plotMarkerSize,marker='o',color='blue',label="Cold Particle")
# Set the axis limits
ax.set_xlim([0,box_size]) 
ax.set_ylim([0,box_size])
ax.axis('equal')
# Set the axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')


if not os.path.exists("movie_frames"):
    os.makedirs("movie_frames")

    
# Function to sequentially load data from each timestep and produce a figure
def update(frame):
    startCold = frame     * nCold
    endCold   = startCold + nCold +1
    startHot  = frame     * nHot
    endHot    = startHot  + nHot+1

    scatterCold.set_offsets(dfCold.iloc[startCold:endCold, :].values)
    scatterHot.set_offsets(dfHot.iloc[startHot:endHot, :].values)
    fig.savefig(f'movie_frames/frame_{frame:04d}.png',dpi=plotDPI)
    return scatterHot,scatterCold



# Create the sequence of figures: (1000 frames takes 11 seconds)
tic = time.time()

#Using parallelization 
threads = 6
with multiprocessing.Pool(processes=threads) as pool:
    pool.map(update,range(timeSteps))    
toc = time.time()

print("Movie frames saved after " + str(round(toc-tic)) + " seconds")
print("\n")
print("\n")

#Generate a movie from the figures using FFmpeg

if not os.path.exists("movie"):
    os.makedirs("movie")

ffmpeg_command = [
    "ffmpeg",
    "-framerate", "300",
    "-i", "movie_frames/frame_%04d.png",
    "-c:v", "libx264",
    "-r", "30",
    "-pix_fmt", "yuv420p",
    "movie/particle-movie.mp4",
    "-t", "1", "-y"
]
subprocess.run(ffmpeg_command)

#All done :)

