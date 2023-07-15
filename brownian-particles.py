import numpy as np

from matplotlib import pyplot as plt
from simfunctions import *
import random
import time
from matplotlib.animation import FuncAnimation
start = time.time()

R = 1e-6
kb = 1.38*10**-23
viscosity = 1*10**-3
#D_T = 1 * 1e-13
D_T = kb*300/(6*np.pi*viscosity*R)
D_R = 1
v = 2e-6
v0 = 50 * 1e-6
delta_t = 0.03
variance = 1
parameters = np.array([R, D_T, D_R, v, delta_t, variance])
n_particles = 100
n_timesteps = 100
box_size = 100 * 1e-6
cut_off_distance = 8 * R
phoreticConstant = v0 * R ** 2
twoPi = 2*np.pi
twoR = 2*R

#Pre-compute certain constants for simulate_motion:
omega = np.sqrt(2*D_T*delta_t)
beta = np.sqrt(2*D_R*delta_t)

position_list = np.zeros([n_particles, 2])
phi_list = np.zeros([n_particles])


def get_distances(_position_candidate, _position_list):
    _distances = np.linalg.norm(_position_candidate - _position_list, axis=1)
    return _distances



def evolve_in_time(_position_list, _phi_list, _n_particles, _box_size, _v0, _cut_off_distance, _n_timesteps,
                   _parameters):
    _R, omega, beta, _v, _delta_t, _variance = parameters
    _position_array = np.zeros([_n_particles, 2, _n_timesteps + 1])
    _discontinuity = np.ones([_n_particles, _n_timesteps + 1])
    _discontinuity_pos = np.ones(_n_particles)
    _position_array[:, :, 0] = _position_list

    for t in range(_n_timesteps):
        if t%10 == True:
          print(f'Current timestep: {t}')
        
        for i in range(_n_particles):
            _phoretic_velocity = np.array([0.0, 0.0])
            for j in range(_n_particles):
                _phoretic_velocity += phoretic_interaction(_v0, _R, _position_array[i, :, t - 1],
                                                           _position_array[j, :, t - 1],
                                                           _cut_off_distance)
            _position_list[i, :] = simulate_motion(_position_list[i, :], _phi_list[i], _phoretic_velocity, _delta_t,
                                                   omega, beta, _v, _v0, _variance, _n_timesteps)

        _position_list, _discontinuity_pos = outside_box(_position_list, _n_particles, _box_size)
        _position_list = check_overlap(_position_list, _R, _n_particles)
        _position_array[:, :, t + 1] = _position_list
        _discontinuity[:, t + 1] = _discontinuity_pos

    return _position_array, _discontinuity


position_list = initialize_positions(position_list, n_particles, box_size)
phi_list = initialize_angles(phi_list, n_particles)
position_list = check_overlap(position_list, R, n_particles)
position_array, discontinuity = evolve_in_time(position_list, phi_list, n_particles, box_size, v0, cut_off_distance, n_timesteps,
                                parameters)

fig1, ax1 = plt.subplots(1, figsize=(8, 8))
disc_vector = np.linspace(0, 2 * np.pi, 100)
#Rescale to fit micrometer scale
position_array *= 1e6
box_size *= 1e6
R *= 1e6

position_array_x = position_array[:, 0, :]
position_array_x = np.multiply(position_array_x, discontinuity)
position_array_y = position_array[:, 1, :]
position_array_y = np.multiply(position_array_y, discontinuity)


for i in range(n_particles):
    #Plot line of particle trajectories
    ax1.plot(position_array_x[i, :], position_array_y[i, :], 'k', linewidth=0.7, alpha=0.1)
    #Plot last particles position
    ax1.plot(position_array_x[i, -1] + R * np.cos(disc_vector), position_array_y[i, -1] + R * np.sin(disc_vector),
             c='k')
    #Plot first particles position
    ax1.plot(position_array_x[i, 0] + R * np.cos(disc_vector), position_array_y[i, 0] + R * np.sin(disc_vector),
             c='k', alpha=0.3)


end = time.time()

print("Execution time = "+str(end-start) + " seconds." )
ax1.set_aspect('equal')
ax1.set_xlabel('x [$\\mu$m]')
ax1.set_xlim([0, box_size])
ax1.set_ylabel('y [$\\mu$m]')
ax1.set_ylim([0, box_size])
plt.tight_layout()
plt.show()

fig, ax = plt.subplots()
ax.set_xlim(0, box_size)  # Set the x-axis limits according to your needs
ax.set_ylim(0, box_size)  # Set the y-axis limits according to your needs
scat = ax.scatter([], [])
#%%
## Movie time:
def update(frame):
    # Update the scatter plot with new positions
    x = position_array_x[:, frame]
    y = position_array_y[:, frame]
    scat.set_offsets(np.column_stack((x, y)))

    return scat,

# Create the animation
animation = FuncAnimation(fig, update, frames=n_timesteps, interval=150)
# Set up the plot window for saving the animation
fig.set_size_inches(8,6)  # Adjust the window size according to your needs
animation.save('particle_animation.gif', writer='pillow')
print("Finished saving .gif")