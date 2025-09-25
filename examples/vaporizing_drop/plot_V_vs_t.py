# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import math

# Load file
data = np.loadtxt('./monitor/simulation', skiprows=2)

# Extract time and volume
t     = data[:, 1 ]
V_num = data[:, 10]
V_ext = data[:, 15]
V_0   = V_ext[0]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Create a figure
fig, ax = plt.subplots(1, 1, figsize=(4, 4))

# Plot
plt.plot(t, V_ext/V_0, ls='-'  , lw=2, color='k')
plt.plot(t, V_num/V_0, ls='--' , lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t (s)$', fontsize=12)
plt.ylabel(r'$V/V_{0}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$Numerical$'], frameon=False, loc='lower left', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./V_vs_t.pdf')