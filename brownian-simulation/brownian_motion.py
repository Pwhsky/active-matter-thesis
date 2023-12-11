import pandas as pd
import time
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
import os
pi	= 3.14159265358979323846
temperature            = 300
particleRadius         = 1e-6 	
cutoff_distance	       = 8*1e-6
delta_t		           = 0.001
kb      		       = 1.380649 * 1e-23
viscosity              = 1e-3
stokesCoefficient	   = 6.0*pi*viscosity*particleRadius
D_T 	               = kb *temperature / stokesCoefficient
trans 	               = np.sqrt(2.0*D_T*delta_t)*0.1

position = [0.0,0.0]


