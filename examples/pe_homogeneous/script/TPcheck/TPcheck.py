

fileName = "Reactor_Log_6-25-2024_new.csv"
newFileName = "Reactor_TPcheck.csv"

f = open(fileName, 'r+')
fnew = open(newFileName, 'w+')

lines = f.readlines()

Time = []
Ts_real = []
Ps_real = []

initTime = ((lines[2].split(',')[0]).split(' ')[-1]).split(':')
initTime = [int(initTime[0]), int(initTime[1]), int(initTime[2])]
initTime = initTime[0]*3600 + initTime[1]*60 + initTime[2]

for i in range(len(lines)):
    if i < 2:
        continue
    else:
        line = lines[i].split(',')
        try:
            Ts_real.append(float(line[-2]))
            Ps_real.append(float(line[-1]))
        except:
            if len(Ts_real) > len(Ps_real):
                Ts_real.pop(-1)
            elif len(Ts_real) < len(Ps_real):
                Ps_real.pop(-1)
            continue

        curTime = ((lines[i].split(',')[0]).split(' ')[-1]).split(':')
        curTime = [int(curTime[0]), int(curTime[1]), int(curTime[2])]
        curTime = curTime[0]*3600 + curTime[1]*60 + curTime[2]

        Time.append(curTime - initTime)

# check length
print(len(Time))
print(len(Ts_real))
print(len(Ps_real))

import CoolProp as CP
import numpy as np
from CoolProp.CoolProp import PropsSI

pi = np.pi

# Define mixture
tankD = 5.72    # cm
tankH = 11.43   # cm

solventVolme = 100.7 # cm^3, initial volume of solvent
NitrogenVolume = pi*tankD**2*tankH/4 - solventVolme # cm^3, initial volume of nitrogen

components = ['nitrogen', 'water']

mole_fraction_liquid = np.zeros(2)
mole_fraction_liquid[0] = 0.0
mole_fraction_liquid[1] = 1.0 - mole_fraction_liquid[0]


def equilibrium(xvals, pressure, components, mole_fraction_liquid):
    '''Function for finding equilibrium temperature and vapor mole fractions.
    
    xvals[0]: temperature (K)
    xvals[1:]: vapor mole fractions
    '''
    temp = xvals[0]
    mole_fraction_gas = [xvals[1], xvals[2]]
    
    pressure_sat = np.zeros(2)
    pressure_sat[0] = PropsSI('P', 'T', temp, 'Q', 1.0, components[0])
    pressure_sat[1] = PropsSI('P', 'T', temp, 'Q', 1.0, components[1])
    
    return [
        (mole_fraction_liquid[0] * pressure_sat[0] - 
         mole_fraction_gas[0] * pressure
         ),
        (mole_fraction_liquid[1] * pressure_sat[1] - 
         mole_fraction_gas[1] * pressure
         ),
        mole_fraction_gas[0] + mole_fraction_gas[1] - 1.0
    ]


print(f'Equilibrium temperature: {sol.x[0]: .2f} K')
print(f'Gas mole fraction of {components[0]}: {sol.x[1]: .3f}')
print(f'Gas mole fraction of {components[1]}: {sol.x[2]: .3f}')


for i in range(len(Time)):
    fnew.write(str(Time[i]) + ',' + str(Ts_real[i]) + ',' + str(Ps_real[i]) + '\n')






f.close()
fnew.close()