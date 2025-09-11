# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Load file
data_nga = np.loadtxt('./monitor/simulation', skiprows=2)
data_tng = np.loadtxt('./R_vs_t_Tanguy')

# Extract data
t_nga = data_nga[:, 1 ]
R_nga = data_nga[:, 11]
R_nga=R_nga/R_nga[0]
R_nga=R_nga**2
t_tng = data_tng[:, 0]
R_tng = data_tng[:, 1]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_tng, R_tng, marker='s', ls='none', color='k')
plt.plot(t_nga, R_nga, ls='-', lw=2, color='b')
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$(R/R_0)^2$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Tanguy ~ et ~ al ~ 2007$', r'$NGA2$'], frameon=False, loc='lower left', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./R_vs_t.pdf')