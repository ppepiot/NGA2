### This project is to simulate the chemical kinetics of polymer pyrolysis

1. convert kinetics file to fortran\
The kinetics file should be from Faravelli group's github repo\
The main script is `script/F.poly2nga2.py`\
usage example:\
`python3 script/F.poly2nga2.py --liquid=PS.liquid --gas=PS.gas --thermo=PS.tdc`\
Note : There are format error in kinetics file, requiring manually corrected.

2. input file\
Number Average Molecular weight (Mn) and Degree of polymerization (DP) decides the ratio between numbers of mid chain roup and end chain group.\
They are from polymer product propertit. \
Temperature can be set as fixed or from temperature.f90 file. The file is created by `script/TPcheck/print_T.py`\
The script requires a temperature log data in .csv format.\
The output interval is set as 100*dt. This can be changed in source code.

3. output file\
The output data are in `monitor`.\
The concentration only includes liquid-phase species. (Assumption: gas-phase species will not participate in reaction)\
In the file mass_fraction, `sum of residure` is the sum of mass change rate of all species. `Liquid mass fraction` is the ratio between mass of all liquid species and initial polymer mass. `LO` stands for light olefin.
