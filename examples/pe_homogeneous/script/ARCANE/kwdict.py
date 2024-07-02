"""
Storing all the possible names for a keyword

"""

import ARCANE.display as display

logger = display.Logger()
logger.set_log('logCompute')


class Kwdict:

    def __init__(self):
        """Class with information on keyword arguments different possible names"""
        # Names (first item is the variable full name)
        self.names = dict()
        self.names['T'] = ['Temperature', 'temperature', 'T']
        self.names['1000/T'] = ['1000/T']
        self.names['HR'] = ['Heat release rate', 'Heat Release', 'Heat release', 'HeatRelease', 'heatrelease',
                            'heat release', 'HR', 'Hr', 'hr',
                            'Heat Release Rate', 'Heat release rate', 'HeatReleaseRate', 'heatreleaserate',
                            'heat release rate', 'HRR', 'Hrr', 'hrr']
        self.names['tig'] = ['Ignition delay time', 'ignition delay time', 'tig', 'ignition']
        self.names['sl'] = ['Laminar flame speed', 'laminar flame speed', 'Sl', 'sl']
        self.names['thickness'] = ['Thickness', 'thickness', 'thick', 'delta', 'd']
        self.names['strain_rate'] = ['Strain Rate', 'strain rate', 'a', 'Strain rate', 'strain_rate']
        self.names['scalar_dissipation_rate'] = ['Scalar Dissipation Rate', 'scalar dissipation rate', 'chi',
                                                 'Scalar dissipation rate', 'scalar_dissipation_rate']
        self.names['Velocity'] = ['Velocity', 'velocity', 'U', 'u']
        self.names['Density'] = ['Density', 'density', 'rho', 'Rho']
        self.names['nu'] = ['Viscosity', 'viscosity', 'nu']
        self.names['Position'] = ['Position', 'position', 'x']
        self.names['Time'] = ['Time', 'time', 't']
        self.names['grid'] = ['Grid', 'grid']
        self.names['grid_interval'] = ['Grid interval', 'grid interval']
        self.names['delta_t'] = ['Delta t', 'delta t', 'dt']
        self.names['delta_x'] = ['Delta x', 'delta x', 'dx']
        self.names['P'] = ['Pressure', 'pressure', 'P', 'p']
        self.names['atom'] = ['Atoms density', 'atom density', 'atoms', 'atom']
        self.names['Y'] = ['Mass fraction', 'MassFraction', 'Mass Fraction', 'massfraction', 'mass fraction', 'Y ', 'Y']
        self.names['X'] = ['Mole fraction', 'MoleFraction', 'Mole Fraction', 'molefraction', 'mole fraction', 'X ', 'X']
        self.names['c'] = ['Concentration', 'concentration', 'c ', 'c']
        self.names['cdot'] = ['Net production rates', 'Production Rates', 'production rates', 'productionrates',
                              'ProductionRates', 'CDOT', 'cdot']
        self.names['w'] = ['Net reaction rates', 'Reaction Rates', 'reaction rates', 'reactionrates',
                           'RRATE', 'WREAC', 'rrate', 'W ', 'w ']
        self.names['Mixture Fraction'] = ['Mixture fraction', 'MixtureFraction', 'Mixture Fraction',
                                          'mixturefraction', 'mixture fraction', 'Z', 'z']
        self.names['phi'] = ['Equivalence ratio', 'Equivalence Ratio', 'equivalence ratio', 'phi', 'Phi']
        self.names['SootVolumeFraction'] = ['Soot volume fraction', 'SootVolumeFraction', 'fv', 'soot_fv']
        self.names['SootNumberDensity'] = ['Soot number density', 'SootNumberDensity', 'Np', 'soot_Np']
        self.names['Xfuel'] = ['Fuel ratio', 'Mole fuel ratio', 'Volume fuel ratio']
        self.names['Yfuel'] = ['Mass fuel ratio']
        self.names['d'] = ['Diameter', 'diameter']

        # Units dictionary
        self.units = dict()
        self.units['Temperature'] = 'K'
        self.units['Position'] = 'm'
        self.units['Diameter'] = 'm'
        self.units['Grid'] = '?'
        self.units['Time'] = 's'
        self.units['Mixture fraction'] = '-'
        self.units['Pressure'] = 'Pa'
        self.units['1000/T'] = '1/K'
        self.units['Equivalence ratio'] = '-'
        self.units['Atoms density'] = '$atom/mol$'
        self.units['Mass fraction'] = '$kg/kg$'
        self.units['Mole fraction'] = '$mol/mol$'
        self.units['Concentration'] = '$kmol/m^3$'
        self.units['Fuel composition'] = '-'
        self.units['Fuel ratio'] = '-'
        self.units['Heat release rate'] = '$J/m^3/s$'
        self.units['Production Rates'] = '$kmol/m^3/s$'
        self.units['Reaction Rates'] = '$kmol/m3/s$'
        self.units['Velocity'] = 'm/s'
        self.units['Viscosity'] = '$m^2/s$'
        self.units['Density'] = '$kg/m^3$'
        self.units['Ignition delay time'] = 's'
        self.units['Laminar flame speed'] = 'm/s'
        self.units['Thickness'] = 'm'
        self.units['integralC0DV'] = 's'
        self.units['integralC0DP'] = 's'
        self.units['integralC1DP'] = 'm'
        self.units['integralC1DCF'] = 'm'
        self.units['Soot volume fraction'] = '-'
        self.units['Soot number density'] = '$cm^{-3}$'
        self.units[' '] = '?'  # Default

        # Reactor label
        self.reactor_labels = dict()

        self.reactor_labels['C0DV'] = ['0D Isochoric reactor', '0DIsochoric', 'C0DV']
        self.reactor_labels['C0DP'] = ['0D Isobaricreactor', '0DIsobaric', 'C0DP']
        self.reactor_labels['C1DP'] = ['1D Premixed freely propagating flame',
                                       '1DPremixed', 'C1DP']
        self.reactor_labels['C1DCFD'] = ['1D Counter flow diffusion flame',
                                         '1DCounterFlowDiffusion', 'C1DCFD']
        self.reactor_labels['C1DF'] = ['1D Flamelet', '1DFlamelet', 'C1DF']
        self.reactor_labels['CXC'] = ['Custom', 'custom', 'CXC']
        self.reactor_labels['C1DB'] = ['1D Burner flame', '1DBurner', 'C1DB']
        self.reactor_labels['C1DCFP'] = ['1D Counter flow premixed flame', '1DCounterFlowPremixed', 'C1DCFP']
        # The following reactors are missing from the new implementation
        self.reactor_labels['C1DIJ'] = ['1D Impinging Jet flame', '1DImpingingJet', 'C1DIJ']
        self.reactor_labels['C1DPSOOT'] = ['1D Premixed freely propagating sooting flame', '1DPremixedSoot', 'C1DPSOOT']
        self.reactor_labels['C1DCFSOOT'] = ['1D Counter flow diffusion sooting flame', '1DCounterFlowSoot', 'C1DCFSOOT']
        self.reactor_labels['C1DBSOOT'] = ['1D Burner sooting flame', '1DBurnerSoot', 'C1DBS']
        self.reactor_labels['C1DIJSOOT'] = ['1D Impinging jet sooting flame', '1DImpingingJetSoot', 'C1DIJSOOT']

        # Reactor grid name
        self.reactor_function_name = dict()
        self.reactor_function_name['C0DV'] = 'Reactor'
        self.reactor_function_name['C0DP'] = 'ConstantPressureReactor'
        self.reactor_function_name['C1DP'] = 'FreeFlame'
        self.reactor_function_name['C1DCFD'] = 'CounterflowDiffusionFlame'
        self.reactor_function_name['C1DF'] = 'Flamelet'
        self.reactor_function_name['CXC'] = ''
        self.reactor_function_name['C1DB'] = 'BurnerFlame'
        self.reactor_function_name['C1DCFD'] = 'CounterflowPremixedFlame'

        # Reactor grid name
        self.reactor_grid_name = dict()
        self.reactor_grid_name['C0DV'] = 'Time'
        self.reactor_grid_name['C0DP'] = 'Time'
        self.reactor_grid_name['C1DP'] = 'Position'
        self.reactor_grid_name['C1DCFD'] = 'Position'
        self.reactor_grid_name['C1DF'] = 'Mixture Fraction'
        self.reactor_grid_name['CXC'] = 'Grid'
        self.reactor_grid_name['C1DB'] = 'Position'
        self.reactor_grid_name['C1DCFP'] = 'Position'
        self.reactor_grid_name['C1DIJ'] = 'Position'
        self.reactor_grid_name['C1DPSOOT'] = 'Position'
        self.reactor_grid_name['C1DCFSOOT'] = 'Position'
        self.reactor_grid_name['C1DBSOOT'] = 'Position'
        self.reactor_grid_name['C1DIJSOOT'] = 'Position'

        # Case arguments needed
        self.required_args = dict()
        self.required_args['C0DV'] = ['pressure', 'temperature']
        self.required_args['C0DP'] = ['pressure', 'temperature']
        self.required_args['C1DP'] = ['pressure', 'temperature']
        self.required_args['C1DCFD'] = ['pressure', 'fuel_temperature']
        self.required_args['C1DF'] = ['pressure', 'fuel_temperature', 'scalar_dissipation_rate']
        self.required_args['C1DB'] = ['pressure', 'temperature']
        self.required_args['C1DIJ'] = ['pressure', 'temperature', 'surface_temperature']
        self.required_args['CXC'] = ['function']
        self.required_args['C1DCFP'] = ['pressure', 'temperature', 'strain_rate']
        self.required_args['C1DPSOOT'] = ['pressure', 'temperature', 'soot_precursors']
        self.required_args['C1DCFSOOT'] = ['pressure', 'fuel_temperature', 'soot_precursors']
        self.required_args['C1DBSOOT'] = ['pressure', 'temperature', 'soot_precursors']
        self.required_args['C1DIJSOOT'] = ['pressure', 'temperature', 'surface_temperature', 'soot_precursors']

        self.optional_args = dict()
        options_base = ['fuel', 'oxidizer', 'phi', 'composition', 'targets', 'error_dict']
        options_1D_specific = ['refine_param', 'ratio', 'slope', 'curve', 'prune', 'width', 'timestep',
                               'imposed_profile']
        options_soot_specific = ['soot_sections', 'soot_precursors']

        self.optional_args['C0DV'] = options_base
        self.optional_args['C0DP'] = options_base
        self.optional_args['C1DP'] = options_base + options_1D_specific + options_soot_specific
        self.optional_args['C1DCFD'] = options_base + options_1D_specific + options_soot_specific \
                                       + ['oxidizer_temperature', 'oxidizer_mass_flow_rate']
        self.optional_args['C1DF'] = options_base + options_1D_specific + options_soot_specific \
                                     + ['oxidizer_temperature', 'oxidizer_mass_flow_rate']
        self.optional_args['C1DB'] = options_base + options_1D_specific + options_soot_specific \
                                     + ['velocity', 'mass_flow_rate']
        self.optional_args['C1DCFP'] = options_base + options_1D_specific + options_soot_specific
        self.optional_args['C1DIJ'] = options_base + options_1D_specific + options_soot_specific \
                                      + ['velocity', 'mass_flow_rate']
        self.optional_args['C1DPSOOT'] = self.optional_args['C1DP']
        self.optional_args['C1DCFSOOT'] = self.optional_args['C1DCFD']
        self.optional_args['C1DBSOOT'] = self.optional_args['C1DB']
        self.optional_args['C1DIJSOOT'] = self.optional_args['C1DIJ']

        self.optional_args['CXC'] = []

        # Reactors default parameters
        self.reactor_default_values = {}

        self.reactor_default_values['C1DP'] = {}
        self.reactor_default_values['C1DP']['width'] = 0.02
        self.reactor_default_values['C1DP']['timestep'] = 1e-9
        self.reactor_default_values['C1DP']['timestep_list'] = [10, 20, 80, 100]

        self.reactor_default_values['C1DCFD'] = {}
        self.reactor_default_values['C1DCFD']['width'] = 0.05
        self.reactor_default_values['C1DCFD']['timestep'] = 1e-9
        self.reactor_default_values['C1DCFD']['timestep_list'] = [1, 2, 5, 10, 20, 80, 100]

        self.reactor_default_values['C1DF'] = {}
        self.reactor_default_values['C1DF']['width'] = 1.0
        self.reactor_default_values['C1DF']['timestep'] = 1e-9
        self.reactor_default_values['C1DF']['timestep_list'] = [10, 20, 80, 100]

        self.reactor_default_values['C1DB'] = {}
        self.reactor_default_values['C1DB']['width'] = 0.02
        self.reactor_default_values['C1DB']['timestep'] = 1e-9
        self.reactor_default_values['C1DB']['timestep_list'] = [1, 2, 5, 10, 20, 80, 100]

        self.reactor_default_values['C1DCFP'] = {}
        self.reactor_default_values['C1DCFP']['width'] = 0.1
        self.reactor_default_values['C1DCFP']['timestep'] = 1e-9
        self.reactor_default_values['C1DCFP']['timestep_list'] = [10, 20, 80, 100]

        # Expendable variables
        self.expandable = ['pressure', 'temperature', 'phi', 'fuel_temperature', 'oxidizer_temperature', 'strain_rate',
                           'fuel_mass_flow_rate', 'oxidizer_mass_flow_rate', 'fuel_velocity', 'oxidizer_velocity',
                           'scalar_dissipation_rate', 'velocity', 'fuel_air_ratio',
                           'mass_flow_rate']  # , 'fuel', 'oxidizer', 'composition']

        # Methods applied on a profile
        self.methods = dict()
        self.methods['max'] = ['Maximum', 'maximum', 'max', 'Max']
        self.methods['min'] = ['Minimum', 'minimum', 'min', 'Min']
        self.methods['mean'] = ['Mean', 'mean', 'Average', 'average']
        self.methods['init'] = ['Initial', 'initial', 'init', 'Init', 'Inlet', 'inlet']
        self.methods['end'] = ['Final', 'final', 'End', 'end']
        self.methods['int'] = ['Integral', 'integral', 'int', 'Int']
        self.methods['dist'] = ['Distance', 'distance', 'dist', 'Dist']
        self.methods['abs'] = ['Absolute', 'absolute', 'abs', 'Abs']

        # Computing methods for scalar variables
        self.compute_methods = {}
        self.compute_methods['sl'] = ['fuel'] + self.methods['abs'] + self.methods['init'] + self.names['HR']
        self.compute_methods['thickness'] = ['thermal', 'diffusive', 'Blint', 'global']
        self.compute_methods['tig'] = self.names['HR']

        # Absolute or relative characterisation
        self.relativity = {}
        self.relativity['error type'] = ['Error type', 'error type', 'Error Type', 'Type', 'type', 'Error', 'error']
        self.relativity['absolute'] = self.methods['abs']
        self.relativity['relative'] = ['Relative', 'relative', 'rel', 'Rel']

        # Name of elements
        # TODO update names an units
        self.elements_names = {}
        # Reducers
        self.elements_names['H'] = 'hydrogen'
        self.elements_names['C'] = 'carbon'
        self.elements_names['Na'] = 'sodium'
        self.elements_names['Al'] = 'aluminium'
        self.elements_names['P'] = 'phosphorus'
        self.elements_names['S'] = 'sulfur'
        self.elements_names['K'] = 'potassium'
        self.elements_names['Fe'] = 'iron'
        self.elements_names['Cu'] = 'copper'
        # Oxidizers
        self.elements_names['O'] = 'oxygen'
        self.elements_names['F'] = 'fluorine'
        self.elements_names['Cl'] = 'chlorine'
        # Inerts
        self.elements_names['N'] = 'nitrogen'
        self.elements_names['Ar'] = 'argon'
        self.elements_names['He'] = 'helium'

        # Weight of elements (kg/mol)
        self.elements_weights = {}
        # Reducers
        self.elements_weights['H'] = 1.00797e-3
        self.elements_weights['C'] = 12.011e-3
        self.elements_weights['Na'] = 22.98977e-3
        self.elements_weights['Al'] = 26.98977e-3
        self.elements_weights['P'] = 30.97376e-3
        self.elements_weights['S'] = 32.06e-3
        self.elements_weights['K'] = 39.97376e-3
        self.elements_weights['Fe'] = 55.845e-3
        self.elements_weights['Cu'] = 63.546e-3
        # Oxidizers
        self.elements_weights['O'] = 15.9994e-3
        self.elements_weights['F'] = 18.998403e-3
        self.elements_weights['Cl'] = 35.453e-3
        # Inerts
        self.elements_weights['N'] = 14.0067e-3
        self.elements_weights['Ar'] = 39.948e-3
        self.elements_weights['He'] = 4.0026e-3

        # Coefficients in the equivalence ratio formula
        # Reducers are taken in a global reaction with oxygen
        # and oxidizers in a global reaction with hydrogen (+ consistency with oxygen)
        self.elements_phi_coeff = {}
        # Reducers
        self.elements_phi_coeff['H'] = 0.5  # H + (1/2) O = (1/2) H2O [1]
        self.elements_phi_coeff['C'] = 2  # C + 2 O = (1) CO2
        self.elements_phi_coeff['Na'] = 1  # Na + O = (1) NaO
        self.elements_phi_coeff['Al'] = 3/2  # Al + 3/2 O = (1/2) Al2O3
        self.elements_phi_coeff['P'] = 5/2  # P + 5/2 O = (1/2) P2O5
        self.elements_phi_coeff['S'] = 2  # S + 2 O = (1) SO2
        self.elements_phi_coeff['K'] = 1  # K + (1) O = (1/2) K2O2
        self.elements_phi_coeff['Fe'] = 3/2  # Fe + 3/2 O = (1/2) Fe2O3
        self.elements_phi_coeff['Cu'] = 1  # Cu +  O = (1) CuO
        # Oxidizers
        self.elements_phi_coeff['O'] = -1  # By definition as the reducer is taken equal to 1
        self.elements_phi_coeff['F'] = -0.5  # H + (1) F = (1) HF with coefficients from [1] 1/2 H + (1/2) F = (1/2) HF
        self.elements_phi_coeff['Cl'] = -0.5  # as for fluorine, 1/2 H + (1/2) Cl = (1/2) HCl
        # Inerts
        self.elements_phi_coeff['N'] = 0
        self.elements_phi_coeff['Ar'] = 0
        self.elements_phi_coeff['He'] = 0


        # # The coefficients can also be determined by the number of valence electrons
        # self.elements_phi_coeff = {}
        # # Reducers
        # self.elements_phi_coeff['H'] = 1
        # self.elements_phi_coeff['C'] = 4
        # self.elements_phi_coeff['Na'] = 1
        # self.elements_phi_coeff['Al'] = 3
        # self.elements_phi_coeff['P'] = 3
        # self.elements_phi_coeff['S'] = 2
        # self.elements_phi_coeff['K'] = 1
        # self.elements_phi_coeff['Fe'] = 2
        # self.elements_phi_coeff['Cu'] = 2
        # # Oxidizers
        # self.elements_phi_coeff['O'] = -2
        # self.elements_phi_coeff['F'] = -1
        # self.elements_phi_coeff['Cl'] = -1
        # # Inerts
        # self.elements_phi_coeff['N'] = 0 * 3
        # self.elements_phi_coeff['Ar'] = 0
        # self.elements_phi_coeff['He'] = 0

    def update(self, attribute_name, new_dictionary):
        """Updates the Kwdict with an external dictionary

        Parameters
        ----------
        attribute_name :
            name of the attribute to update
        new_dictionary :
            dictionary to be added

        """

        if hasattr(self, attribute_name):
            temp_dict = getattr(self, attribute_name)
            temp_dict.update(new_dictionary)

            setattr(self, attribute_name, temp_dict)
        else:
            logger.error("The Keyword dictionary cannot be updated has it has no attribute: " + attribute_name)

    def get_full_name(self, name):
        """Get the full name of abbreviated parameter

        Parameters
        ----------
        name :
            abbreviated name of the parameter

        Returns
        -------
        full_name : str
            full name of the parameter

        """

        full_name = ''
        for key in self.names:
            if name in self.names[key]:
                full_name = self.names[key][0]

        for key in self.methods:
            if name in self.methods[key]:
                full_name = self.methods[key][0]

        if not full_name:
            logger.debug('No match found for full name')
            full_name = ' '

        return full_name

    def get_base_name(self, name):
        """Get the base name of abbreviated parameter

        Parameters
        ----------
        name :
            abbreviated name of the parameter

        Returns
        -------
        base_name : str
            base name of the parameter

        """

        base_name = ''
        for key in self.names:
            if name in self.names[key]:
                base_name = key
                break

        if not base_name:
            logger.debug('No match found for base name')
            base_name = name

        return base_name

    def get_reactor_name(self, name, short=False):
        """Get the full name of a reactor from the abbreviated name

        Parameters
        ----------
        name :
            abbreviated name
        short :
            if True, gets the short name (Default value = False)

        Returns
        -------
        full_name : str
            full name of the reactor

        """

        full_name = ''
        if short:
            name_index = -1
        else:
            name_index = 0

        for key in self.reactor_labels:
            if name in self.reactor_labels[key]:
                full_name = self.reactor_labels[key][name_index]

        if not full_name:
            logger.debug('No match found for full name')
            full_name = ' '

        return full_name

    def expand_names(self, data_dict):
        """Expands the names dictionary with kwdict names entries

        Parameters
        ----------
        data_dict :
            Dictionary containing the data

        Returns
        -------
        data_dict : dict
            New dictionary with expanded entries

        """

        new_data_dict = {}

        for name in data_dict:
            # print('name', name)
            for base_name in self.names:
                # print('    base name', base_name)
                if name in self.names[base_name]:
                    # print('        match', self.names[base_name])
                    for other_name in self.names[base_name]:
                        new_data_dict[other_name] = data_dict[name]
            # input()
        # print(new_data_dict.keys())
        data_dict.update(new_data_dict)

        return data_dict
