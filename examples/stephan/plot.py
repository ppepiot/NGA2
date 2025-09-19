# Import libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import math

# Load file
data = np.loadtxt('./monitor/simulation', skiprows=2)

# Extract data
t = data[:, 1 ]
x = data[:, 11] * 1000
u = data[:, 12]
x_ext = data[:, 13] * 1000
u_ext = data[:, 14]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot
fig, ax = plt.subplots(1, 1, figsize=(4, 4))

plt.plot(t, 1000*x_ext, ls='-', lw=2, color='k')
plt.plot(t, 1000*x    , ls='-', lw=2, color='b')
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$x_{I}~(mm)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$Numerical$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./x_vs_t.pdf')

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t, 1000*u_ext, ls='-', lw=2, color='k')
plt.plot(t[1:], 1000*u[1:], ls='-', lw=2, color='b')
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$u_I~(mm/s)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$Numerical$'], frameon=False, loc='upper right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./u_vs_t.pdf')