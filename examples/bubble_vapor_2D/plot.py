# Import libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import math

# Constants
tI = 0.15216

# Load file
data = np.loadtxt('./monitor/simulation', skiprows=2)

# Extract data
t = data[:, 1 ]
t = tI + t
R = data[:, 11]
R_ext = data[:, 12]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t, R_ext, ls='-', lw=2, color='k')
plt.plot(t, R    , ls='-', lw=2, color='b')
# plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$R~(m)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$Numerical$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./R_vs_t.pdf')