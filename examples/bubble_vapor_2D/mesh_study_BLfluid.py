# Import libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import math

# Constants
tI = 0.5

# Load file
data_128 = np.loadtxt('./monitor/simulation_128_BLfluid', skiprows=2)
data_256 = np.loadtxt('./monitor/simulation_256_BLfluid', skiprows=2)


# Extract data
R_ext  = data_256[:, 12]

t_128  = data_128[:, 1 ]
t_128  = tI + t_128
R_128  = data_128[:, 11]

t_256  = data_256[:, 1 ]
t_256  = tI + t_256
R_256  = data_256[:, 11]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_256, R_ext, ls='-' , lw=2, color='k')
plt.plot(t_128, R_128, ls='-' , lw=2, color='g')
plt.plot(t_256, R_256, ls='-' , lw=2, color='b')
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$R~(m)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$NGA 2 ~ 128 \times 128$', r'$NGA 2 ~ 256 \times 256$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./R_vs_t_BLfluid.pdf')