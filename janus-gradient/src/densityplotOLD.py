import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time 
import subprocess
import sys
#TODO: MAKE A HEATMAP DENSITY PLOT THAT LOOKS SICK AF
df = pd.read_csv('deposits.csv')




num_bins = 50
max_distance = 5e-6
rdf_values = np.zeros(num_bins)

distances = np.linalg.norm(df[['x', 'y', 'z']].values[:, np.newaxis, :] - df[['x', 'y', 'z']].values, axis=2)
np.fill_diagonal(distances, np.inf)

bins = np.linspace(0, max_distance, num_bins + 1)

# Calculate the RDF by counting distances within each bin
for i in range(num_bins):
    bin_mask = (distances >= bins[i]) & (distances < bins[i + 1])
    rdf_values[i] = (bin_mask.sum() / (4 * np.pi / 3 * (bins[i + 1]**3 - bins[i]**3))) / len(df)

# Normalize the RDF by the average density


average_density = len(df) / (4 * np.pi / 3 * 2e-6**3)


rdf_values /= average_density



plt.figure(figsize=(8, 6))
plt.plot(bins[2:-1], rdf_values[1:-1])
plt.plot()
plt.xlabel('Radial distance [m]')
plt.ylabel('RDF')
plt.title('Radial Distribution Function of Deposits')
plt.grid(True)
plt.savefig("density.jpg")
