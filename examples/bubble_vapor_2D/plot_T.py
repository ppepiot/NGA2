# Import libraries
import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Date file reader
def file_to_data(file, data_names):

    # Read lines
    with open(file, 'r') as sf:
        lines = sf.readlines()

    # Store header
    names = lines[0].split()

    # Remove header line
    lines.pop(0)

    fdata = []
    for line in lines:
        split_line = line.split()
        fdata.append([float(val) for val in split_line])
    fdata = np.array(fdata)

    data = {}
    for name in data_names:
        data[name] = fdata[:, names.index(name)]
    
    return data

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Path where the files are located
data_dir = '/Users/shayanhbi/Repositories/NGA2/examples/bubble_vapor_2D/temperature'

# Get and sort data
arr_names = np.array(['x', 'T'])
data_ext = []
data_num = []
time_num = []
time_ext = []

pattern = r'[-+]?\d*\.\d+|\d+'
for filename in os.listdir(data_dir):
    if filename.endswith('.dat'):
        match = re.search(pattern, filename)
        time = float(match.group())
        if 'num' in filename:
            time_num.append(time)
            data_num.append(file_to_data(data_dir +'/' + filename, arr_names))
        elif 'ext' in filename:
            time_ext.append(time)
            data_ext.append(file_to_data(data_dir + '/' + filename, arr_names))
        else:
            pass

sorted_lists = sorted(zip(time_num, data_num))
time_num, data_num = zip(*sorted_lists)

sorted_lists = sorted(zip(time_ext, data_ext))
time_ext, data_ext = zip(*sorted_lists)


# Make colors
colors = plt.cm.plasma(np.linspace(0, 1, len(time_ext)))

fig, ax = plt.subplots(1, 1, figsize=(8,6))

# Visualize
for i, t in enumerate(time_ext):
    ax.plot(data_ext[i]['x'], data_ext[i]['T'], '-',  color=colors[i], label = r'$t = {:.3f}~(s)$'.format(t))
    ax.plot(data_num[i]['x'], data_num[i]['T'], '--', color=colors[i])

# Custom legend
custom_lines = [
    Line2D([0], [0], color='k', lw=1.5, ls='-',  label=r'$Exact$'),
    Line2D([0], [0], color='k', lw=1.5, ls='--', label=r'$Numerical$'),
]
first_legend = ax.legend(custom_lines, [r'$Exact$', r'$Numerical$'], frameon=False, loc='lower right', bbox_to_anchor=(1, 0.3), fontsize=14)
ax.add_artist(first_legend)
ax.legend(frameon=False, loc='lower right', fontsize=14)
plt.grid(which='major', axis='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.25)
plt.xlabel(r'$x$', fontsize=14)
plt.ylabel(r'$T$', fontsize=14)
for spine in plt.gca().spines.values():
    spine.set_linewidth(1.2)
plt.tight_layout()
# plt.xlim(1, 0.4)
# plt.ylim(1, 3)
plt.savefig('./T.pdf')
