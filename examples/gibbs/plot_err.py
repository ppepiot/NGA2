# Import libraries
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('./newton', skiprows=1)
iter = data[:, 0]
err  = data[:, 1]
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(iter,err, ls='-', lw=2, color='k')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# ax.set_yscale('log')
plt.xlabel(r'$Iteration ~ number$', fontsize=12)
plt.ylabel(r'$Residual ~ norm$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
plt.tight_layout()
plt.savefig('./err.pdf')