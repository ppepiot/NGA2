# Import libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import math

# Load file
data = np.loadtxt('./RH_vs_T', skiprows=1)

# Extract data
T  = data[:, 0]
RH = data[:, 1]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(T, RH, ls='-' , lw=2, color='b')
# plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$T$', fontsize=12)
plt.ylabel(r'$R_H$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
plt.tight_layout()
plt.savefig('./RH.pdf')