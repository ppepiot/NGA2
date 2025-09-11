# Import libraries
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Path where your .dat files are located
data_dir = "/Users/shayanhbi/Repositories/NGA2/examples/Stephan1D/temperature"

# Get file lists
ext_files = sorted(glob.glob(os.path.join(data_dir, "ext*.dat")))
nga_files = sorted(glob.glob(os.path.join(data_dir, "NGA2*.dat")))

def get_timestamp(fname, prefix):
    base = os.path.basename(fname)
    return base.replace(prefix, "").replace(".dat", "")

plt.figure(figsize=(8,6))

for extf, ngaf in zip(ext_files, nga_files):
    t = get_timestamp(os.path.basename(extf), "ext")

    ext_data = np.loadtxt(extf, comments='#')
    nga_data = np.loadtxt(ngaf, comments='#')

    plt.plot(1000*ext_data[:, 0], ext_data[:, 1], '-' , color='k')
    plt.plot(1000*nga_data[:, 0], nga_data[:, 1], '--', color='k')

plt.grid(which='major', axis='both', color='gray', linestyle='-', linewidth=0.4, alpha=0.25)
plt.xlabel(r'$x~(mm)$', fontsize=12)
plt.ylabel(r'$T_g~(K)$', fontsize=12)
plt.tight_layout()
plt.savefig('./T_g.pdf')
