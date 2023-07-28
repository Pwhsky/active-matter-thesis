import pandas as pd
import time
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from concurrent.futures import ThreadPoolExecutor
from matplotlib.animation import PillowWriter
from functools import partial
from PIL import Image

#%%

plotMarkerSize = 30
plotDPI        = 200


n_hot       = 20
n_cold      = 30
nParticles = n_hot + n_cold
timeSteps  = 500
#Handler to simulate in C++:
subprocess.run(['g++','brownian-particles.cpp','-o','sim','-O4'])
subprocess.run(['./sim',f' {n_hot}',f' {n_cold}',f' {timeSteps}'])

print("Processing data...")
box_size = 100 * 1e-6

# Read the CSV file into a DataFrame
with open('hotParticles.csv', 'r') as file:
    df_hot = pd.read_csv(file, header=None, names=['x', 'y'],low_memory=True)

with open('coldParticles.csv', 'r') as file:
    df_cold = pd.read_csv(file, header=None, names=['x', 'y'], low_memory=True)

# Extract the x and y coordinates from the DataFrame
df_hot['x'] = pd.to_numeric(df_hot['x'],errors='coerce')
df_hot['y'] = pd.to_numeric(df_hot['y'],errors='coerce')

df_cold['x'] = pd.to_numeric(df_cold['x'],errors='coerce')
df_cold['y'] = pd.to_numeric(df_cold['y'],errors='coerce')



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
def update(frame, df_hot, df_cold, n_hot, n_cold):
    if frame % 3 == 0 and frame % 3 != 1:
        start_cold = frame * n_cold
        end_cold = start_cold + n_cold + 1
        start_hot = frame * n_hot
        end_hot = start_hot + n_hot + 1

        scatterCold.set_offsets(df_cold.iloc[start_cold:end_cold, :].values)
        scatterHot.set_offsets(df_hot.iloc[start_hot:end_hot, :].values)
    return scatterHot, scatterCold

# Function to generate animation frames for a given chunk of frames
def generate_animation_frames(chunk):
    frames = range(chunk * n_hot, (chunk + 1) * n_hot + 1)
    with ThreadPoolExecutor() as executor:
        update_partial = partial(update, df_hot=df_hot, df_cold=df_cold, n_hot=n_hot, n_cold=n_cold)
        results = list(executor.map(update_partial, frames))
    return results

# Number of chunks for parallelization (adjust this based on your CPU cores)
num_chunks = 4

# Create the animation frames in parallel
animation_frames = []
with ThreadPoolExecutor() as executor:
    animation_frames = list(executor.map(generate_animation_frames, range(num_chunks)))

# Flatten the frames list
animation_frames = [frame for frames_chunk in animation_frames for frame in frames_chunk]

# Create the animation
ani = animation.ArtistAnimation(fig, animation_frames, interval=4,blit=True)

# Save each frame as an image
image_frames = []
for frame in ani.new_frame_seq():
    fig.canvas.draw()
    image = Image.frombytes("RGB", fig.canvas.get_width_height(), fig.canvas.tostring_rgb())
    image_frames.append(image)

# Save the animation as an animated GIF
gif_file = "particle_movie.gif"
image_frames[0].save(gif_file, format="GIF", append_images=image_frames[1:], save_all=True, duration=100, loop=0)

# ... (rest of your code)
# Show the plot (optional)
gif_file = "particle_movie.gif"
toc = time.time()



subprocess.run(["xdg-open",gif_file])

#Close and clear stuff to prevent leaks
plt.close(fig)
del ani


