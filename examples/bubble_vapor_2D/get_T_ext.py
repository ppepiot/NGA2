# Import libraries
import os
import math
import numpy as np
from scipy.integrate import quad

def integrand(zeta, beta, rho_g, rho_l):
    return np.exp(-beta**2 * ((1 - zeta)**(-2) - 2 * (1 - rho_g / rho_l) * zeta - 1)
)

def get_integral(beta, rho_g, rho_l, R, r):
    lower_limit = 1 - R / r
    upper_limit = 1.0
    result, err = quad(integrand, lower_limit, upper_limit, args=(beta, rho_g, rho_l))
    return result, err

res_path = '/Users/shayanhbi/Repositories/NGA2/examples/bubble_vapor_2D/temperature'
if not os.path.exists(res_path):
    try:
        os.mkdir(res_path)
    except OSError as error:
        print(f"Error: {error}")

rhog = 0.25
kg = 0.007
Cpg = 1
rhol = 2.5
kl = 0.07
Cpl = 2.1
hlg=100
alphal = kl / (rhol * Cpl)
Tinf = 3
Tsat = 1
beta = 0.69252756051719189
cnst=2*beta**2*(rhog*(hlg+(Cpl-Cpg)*(Tinf-Tsat)))/(rhol*Cpl)

desired_times = [0.5, 0.6, 0.7, 0.8]

nx = 200
x  = np.linspace(0, 0.565, nx)
for t_ind, t in enumerate(desired_times):
    print(t)
    x_i = 2 * beta * np.sqrt(alphal * (t))
    T  = np.zeros(len(x))
    for i in range(len(x)):
        if (x[i] <= x_i):
            T[i] = Tsat
        else:
            int, err = get_integral(beta, rhog, rhol, x_i, x[i])
            T[i] = Tinf - cnst * int
    with open(res_path + '/ext' + str(desired_times[t_ind]) + '.dat', 'w') as file_ctr:
        file_ctr.write(f"x            T\n")
        for i in range(len(x)):
            file_ctr.write(f"{x[i]}            {T[i]}\n")
    file_ctr.close()
