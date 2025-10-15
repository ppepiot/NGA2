# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import math

# Load file
data = np.loadtxt('./monitor/mflux', skiprows=2)

# Extract time and volume
n  = data[:, 0]
eL = data[:, 2]
eG = data[:, 3]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Create a figure
fig, ax = plt.subplots(1, 1, figsize=(4, 4))

# Plot
plt.plot(n, eL, ls='-' , lw=2, color='r')
plt.plot(n, eG, ls='-' , lw=2, color='b')
# plt.yscale('log')
plt.xlabel(r'$Pseudo\,time\,step\,number$', fontsize=12)
plt.ylabel(r'$Mass\,flux\,error$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Liquid$', r'$Gas$'], frameon=False, loc='upper right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./mfllux_err.pdf')