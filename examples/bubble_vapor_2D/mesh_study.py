# Import libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import math

# Constants
tI = 0.15216

# Load file
data_ext = np.loadtxt('./monitor/simulation_ext', skiprows=2)
data_128 = np.loadtxt('./monitor/simulation_128', skiprows=2)
data_256 = np.loadtxt('./monitor/simulation_256', skiprows=2)
data_512 = np.loadtxt('./monitor/simulation_512', skiprows=2)
data_cip = np.loadtxt('./Cipriano_256', skiprows=2)
data_mf  = np.loadtxt('./massFlux_vs_grid')


# Extract data
t_ext  = data_ext[:, 1 ]
R_ext  = data_ext[:, 12]
t_128  = data_128[:, 1 ]
t_128  = tI + t_128
R_128  = data_128[:, 11]
t_256  = data_256[:, 1 ]
t_256  = tI + t_256
R_256  = data_256[:, 11]
t_512  = data_256[:, 1 ]
t_512  = tI + t_512
R_512  = data_256[:, 11]
t_cip  = data_cip[:, 0 ]
R_cip  = data_cip[:, 1 ]
n      = data_mf [:, 0 ]
mf_num = data_mf [:, 1 ]
mf_ext = data_mf [:, 2 ]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_ext, R_ext*1000, ls='-' , lw=2, color='k')
plt.plot(t_128, R_128*1000, ls='-' , lw=2, color='g')
plt.plot(t_256, R_256*1000, ls='-' , lw=2, color='b')
plt.plot(t_512, R_512*1000, ls='--', lw=2, color='r')
plt.plot(t_cip, R_cip*1000, ls='-.', lw=2, color='c')
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$R~(mm)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$NGA 2 ~ 128 \times 128$', r'$NGA 2 ~ 256 \times 256$', r'$NGA 2 ~ 512 \times 512$', r'$Cipriano ~ et ~ al ~ 256 \times 256$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./R_vs_t.pdf')

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(n, mf_ext, ls='-' , lw=2, color='k')
plt.plot(n, mf_num, marker='s', ls='none', color='k')
plt.xlabel(r'$No. ~ grid ~ cells$', fontsize=12)
plt.ylabel(r'$\dot{m}^{\prime\prime}$', fontsize=12)
# ax.set_xticks(n)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
plt.tight_layout()
plt.savefig('./mf_vs_grid.pdf')