# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Load file
dataNGA = np.loadtxt('./T0')
dataBL  = np.loadtxt('./T0_BL')

# Extract data
x  = dataNGA[:, 0]
T0 = dataNGA[:, 1]

x_BL  = dataBL[:, 0]
T0_BL = dataBL[:, 1]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(x, T0, ls='-', lw=2, color='k')
plt.plot(x_BL, T0_BL, marker='s', ls='none', color='b')
plt.xlabel(r'$x~(m)$', fontsize=12)
plt.ylabel(r'$T_0~(K)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$NGA2$', r'$Boyd ~ and ~ Ling$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./T0.pdf')