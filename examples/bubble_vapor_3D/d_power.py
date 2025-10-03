# Import libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import math

# Load file
data_ext     = np.loadtxt('./monitor/simulation_ext', skiprows=2)
data_256_d2  = np.loadtxt('./monitor/simulation_256_dpow2', skiprows=2)
data_256_d_2 = np.loadtxt('./monitor/simulation_256_dpow-2', skiprows=2)
data_256_d_4 = np.loadtxt('./monitor/simulation_256_dpow-4', skiprows=2)


# Extract data
t_ext  = data_ext[:, 1 ]
R_ext  = data_ext[:, 12]

t_256_d2  = data_256_d2[:, 1 ]
R_256_d2  = data_256_d2[:, 11]

t_256_d_2  = data_256_d_2[:, 1 ]
R_256_d_2  = data_256_d_2[:, 11]

t_256_d_4  = data_256_d_4[:, 1 ]
R_256_d_4  = data_256_d_4[:, 11]


# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
# plt.plot(t_ext, R_ext*1000, ls='-' , lw=2, color='k')
plt.plot(t_256_d2, R_256_d2*1000, ls='-' , lw=2, color='b')
plt.plot(t_256_d_2, R_256_d_2*1000, ls='-' , lw=2, color='r')
plt.plot(t_256_d_4, R_256_d_4*1000, ls='-' , lw=2, color='g')
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$R~(mm)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
# legends = ax.legend([r'$Analytical$', r'$NGA 2 ~ 256 \times 256 ~ d^{2}$', r'$NGA 2 ~ 256 \times 256 ~ d^{-2}$', r'$NGA 2 ~ 256 \times 256 ~ d^{-4}$'], frameon=False, loc='lower right', fontsize=10)
legends = ax.legend([r'$NGA 2 ~ 256 \times 256 ~ d^{2}$', r'$NGA 2 ~ 256 \times 256 ~ d^{-2}$', r'$NGA 2 ~ 256 \times 256 ~ d^{-4}$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./d_power.pdf')