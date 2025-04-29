# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import math

# Constants
Lx = 0.001
Ly = 0.001
Lz = 7.8125e-06
H = 1e-4
tI      = 0.027
beta    = 0.067819663373413
alpha_g = 0.025 / (0.597 * 2030)
Tw      = 383.15
Tsat    = 373.15

# Load file
data = np.loadtxt('./monitor/simulation', skiprows=2)

# Extract data
t = data[:, 1 ]
t = tI + t
x = data[:, 11] * 1000
vol_l = data[:, 10]
u = data[:, 4]
print(min(u))

# Analytical solution
x_ext = 2 * beta * np.sqrt(alpha_g * t) * 1000
vol_l_ext = (Lx - H) * Ly * Lz
u_ext = beta * np.sqrt(alpha_g / t)

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot
fig, ax = plt.subplots(1, 1, figsize=(4, 4))

plt.plot(t, x_ext, ls='-', lw=2, color='k')
plt.plot(t, x    , ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t$', fontsize=12)
plt.ylabel(r'$x_{I}~(mm)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$Numerical$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./x_vs_t.pdf')

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t, u_ext, ls='-', lw=2, color='k')
plt.plot(t, u    , ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t$', fontsize=12)
plt.ylabel(r'$u~(m/s)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$Numerical$'], frameon=False, loc='upper right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./u_vs_t.pdf')