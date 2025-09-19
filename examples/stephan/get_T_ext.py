# Import libraries
import os
import math
import numpy as np

res_path = '/Users/shayanhbi/Repositories/NGA2/examples/stephan/temperature'
if not os.path.exists(res_path):
    try:
        os.mkdir(res_path)
    except OSError as error:
        print(f"Error: {error}")

rhog = 0.597
kg = 0.025
Cpg = 2030
alphag = kg / (rhog * Cpg)
Tw = 383.15
Tsat = 373.15
beta = 6.6916063714766549E-002

desired_times = [0.027, 0.1, 0.2, 0.3, 0.4, 0.5]

nx = 100
for t_ind, t in enumerate(desired_times):
    x_i = 2 * beta * np.sqrt(alphag * (t))
    x   = np.linspace(0, x_i, nx)
    Tg  = np.zeros(len(x))
    for i in range(len(x)):
        Tg[i]  = Tw + (Tsat - Tw) / math.erf(beta) * math.erf(x[i] / (2 * np.sqrt(alphag*(t))))
    with open(res_path + '/ext' + str(desired_times[t_ind]) + '.dat', 'w') as file_ctr:
        file_ctr.write(f"x            Tg\n")
        for i in range(len(x)):
            file_ctr.write(f"{x[i]}            {Tg[i]}\n")
    file_ctr.close()