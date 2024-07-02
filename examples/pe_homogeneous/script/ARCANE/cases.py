""" Management of 0D/1D run cases

This module creates and manages all runs necessary for reduction and validation
It includes class definitions for the various types of simulations
* 0D, isochoric or isobaric
* 1D, freely propagating premixed flames

It also includes various functions that act on cases
* etc.

"""
import ARCANE.tools as tools
import ARCANE.kwdict as kwdict
import ARCANE.error as error
import ARCANE.display as display
import ARCANE.postproc as postproc
import ARCANE.database as database

import cantera as ct

import numpy as np
import sys
import csv

from copy import deepcopy

import subprocess
import os.path
import time

logger = display.Logger()
logger.set_log('logCompute')

kwdict = kwdict.Kwdict()

# Avoid thermodynamics warning being printed thousands of times
ct.suppress_thermo_warnings()


def init_case_database(casedir='cases', overwrite=False):
    """Initialize case folder.
    For now, do nothing if it exists. Might include cleanup routines later on.

    Parameters
    ----------
    casedir :
        name of folder that will contain all case simulations, default is 'cases'
    overwrite :
        optional Boolean. if True, folder is overwritten (Default value = False)

 
        
    Created: 17/11/12 [PP]
    
    Last modified: 18/10/02 [QC]

    """

    # Initialize case directory
    database.create_dir(casedir, overwrite)

    # Set up casedir for all instances
    Case.casedir = casedir


def example_case(reactor_type='C0DV', mechanism=None):
    """Generates an example class :func:`~ARCANE.mechanisms.Mechanism` object at classical conditions.

    Parameters
    ----------
    reactor_type :
        type of reactor to be computed (Default value = 'C0DV')
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)

    Returns
    -------
    example class :func:`~ARCANE.mechanisms.Mechanism` object

    """

    if not mechanism:
        import ARCANE.mechanisms as mechanisms
        mechanism = mechanisms.example_mechanism()

    # Changing the temperature between 0D and 1D
    if reactor_type in kwdict.reactor_labels['C0DV'] + kwdict.reactor_labels['C0DP']:
        temperature = 1000
    else:
        temperature = 300

    # Using standard values
    pressure = 1e5
    phi = 1

    # Methane and air are chosen
    fuel = 'X/CH4/1'
    oxidizer = 'X/O2/21/N2/79'

    new_case = create_case(reactor=reactor_type,
                           mechanism=mechanism,
                           fuel=fuel,
                           oxidizer=oxidizer,
                           pressure=pressure,
                           temperature=temperature,
                           phi=phi)[0]

    new_case.run()

    return new_case


def run_cases(cases_list, mechanisms_list, timing=False, overwrite=False, restore='auto',
              permissive=True, check_validity=True):
    """
    Go through all case instances and run case using specified mechanism.

    Parameters
    ----------
    cases_list :
        list of the cases to run
    mechanisms_list :
        list of class :func:`~ARCANE.mechanisms.Mechanism` object containing all info about mechanism
    timing :
        if True, return the time of the simulation as an array of cases_list length (Default value = False)
    overwrite :
        if True, overwrites existing data,
        this argument is directly passed to :func:`~ARCANE.cases.Case.run` method (Default value = False)
    restore :
        if auto, restores the solution of the parent mechanism or closest case,
        this argument is directly passed to :func:`~ARCANE.cases.Case.run` method (Default value = 'auto')
    permissive :
        if False, run_cases will stop on the first failed case (Default value = True)
    check_validity : bool
        if True, checks every case to ensure its success else only computes the not already computed cases


    Created: 17/11/13 [PP]

    Last modified: 19/09/10 [QC]

    """

    loglevel = display.cantera_solve_loglevel

    if not isinstance(cases_list, list):
        cases_list = [cases_list]

    if not isinstance(mechanisms_list, list):
        mechanisms_list = [mechanisms_list]

    sim_time = np.zeros(len(cases_list) * len(mechanisms_list))

    count = 0

    for index_mech, mechanism in enumerate(mechanisms_list):

        if mechanism.f90 and not mechanism.dynamic_lib_opened:
            mechanism.reset()
            mechanism.compile()

        for case in cases_list:
            case.success = False

        if not check_validity:
            # Retrieving the cases already computed
            cases_list_ids = [case.myid for case in cases_list]
            computed_cases_ids = mechanism.computed_cases(create_case=False)
            not_computed_cases = [case for case in cases_list if case.myid not in computed_cases_ids]
            not_computed_cases_ids = [case.myid for case in not_computed_cases]
            reduced_cases_list = not_computed_cases

            logger.info(f"{len(not_computed_cases_ids)} cases need to be computed (out of {len(cases_list_ids)}).")

        else:
            reduced_cases_list = cases_list

        for index_case, case in enumerate(reduced_cases_list):
            start_time = time.time()
            case.run(mechanism, loglevel=loglevel, overwrite=overwrite, restore=restore, permissive=permissive)
            end_time = time.time()

            sim_time[count] = end_time - start_time

            count += 1

            if not case.success and not permissive:
                logger.warning(case.myid + " failed so the running stopped ! Keep the 'permissive' keyword as True "
                                           "to continue the other cases computation")
                break

        if permissive:
            # Checking the success of the cases
            for i, case in enumerate(reduced_cases_list):
                if case.ndata < 0:
                    case.data_dict(mechanism)

                if not case.success or case.ndata < 10:
                    logger.info(f"Case {case.solution} failed on first try. Trying again !")

                    case.run(mechanism, loglevel=loglevel, overwrite=True, restore=False, permissive=permissive)

        if mechanism.f90 and mechanism.dynamic_lib_opened:
            mechanism.reset()

    if not timing:
        sim_time = None

    return sim_time


def create_case(reactor, mechanism, **kwargs):
    """Function creating the case object given input parameters.

    Parameters
    ----------
    reactor :
        type of reactor to be created
    mechanism :
        reference class :func:`~ARCANE.mechanisms.Mechanism` object
    kwargs :
        any argument that can parametrize the case
        

    Returns
    -------
    list of cases generated

    """

    parameters_dict = kwargs
    expandable_keys = [key for key in list(kwargs.keys()) if key in kwdict.expandable]

    fixed_keys = [key for key in list(kwargs.keys()) if key not in kwdict.expandable]
    fixed_dict = dict()
    for key in fixed_keys:
        fixed_dict.update({key: parameters_dict[key]})

    count = 0
    storing_list = list()

    # Expanding values
    for parameter in parameters_dict:
        if parameter in expandable_keys:
            # Converting floats and list to correct format
            if type(parameters_dict[parameter]) in [int, float]:
                parameters_dict[parameter] = [parameters_dict[parameter]]
            elif type(parameters_dict[parameter]) == list:
                parameters_dict[parameter] = parameters_dict[parameter]
            else:
                parameters_dict[parameter] = tools.expand_values(parameters_dict[parameter])

    while count < len(expandable_keys):

        size = len(storing_list)
        new_list = list()

        if size > 0:
            for step in range(size):
                for value in parameters_dict[expandable_keys[count]]:
                    local_dict = storing_list[step].copy()
                    local_dict.update({expandable_keys[count]: value})

                    new_list.append(local_dict)

        else:
            for value in parameters_dict[expandable_keys[count]]:
                new_list.append({expandable_keys[count]: value})

        storing_list = new_list.copy()
        count += 1

    # Creating all cases
    new_cases = list()

    for parameters in storing_list:
        parameters.update(fixed_dict)
        new_cases.append(Case(reactor_type=reactor, mechanism=mechanism, **parameters))

    # Case with not expandable parameters
    if not storing_list:
        new_cases.append(Case(reactor_type=reactor, mechanism=mechanism, **fixed_dict))

    return new_cases


class Case(object):
   

    # Generator
    def __init__(self, reactor_type, mechanism, **kwargs):
        """Main Case class, containing all class-nonspecific methods and attributes.
    
        Constructor is specified based on first attribute being passed.
        
        Constructor of super class Case.
        
        :param reactor: type of reactor to be created
        :param mechanism: reference class :func:`~ARCANE.mechanisms.Mechanism` object
        :param kwargs: any argument that can parametrize the case

        """

        self.myid = None
        self.nickname = None
        self.success = True

        self.ndata = -1
        self.reduced_data = None

        # Equilibrium values
        self.equilibrium = None
        self.temperature_eq = None
        self.pressure_eq = None
        self.composition_eq = None

        # Type of reactor
        found = False
        for reactor_name in kwdict.reactor_labels:
            if reactor_type in kwdict.reactor_labels[reactor_name]:
                reactor_type = reactor_name
                found = True

        possible_names = [y for x in kwdict.reactor_labels.values() for y in x]

        if not found:
            logger.warning('The specified reactor type ' + reactor_type +
                           ' does not exist')
            logger.warning('Available reactor types are: ' + ', '.join(possible_names))
            quit()

        self.reactor_type = reactor_type

        # Reference mechanism
        self.mechanism = mechanism

        # Updating kwdict with the species_names
        kwdict.update('names', self.mechanism.kwdict_names_extension)
        kwdict.update('units', self.mechanism.kwdict_units_extension)

        # Dictionary of keywords
        self.parameters = kwargs

        missing_parameters = [key for key in kwdict.required_args[self.reactor_type]
                              if key not in list(self.parameters.keys())]

        if 'phi' in missing_parameters and ('composition' in list(self.parameters.keys())
                                            or 'fuel_air_ratio' in list(self.parameters.keys())):
            missing_parameters.remove('phi')
            if 'fuel' in missing_parameters:
                missing_parameters.remove('fuel')
            if 'oxidizer' in missing_parameters:
                missing_parameters.remove('oxidizer')

        if 'strain_rate' in missing_parameters \
                and ('fuel_mass_flow_rate' in list(self.parameters.keys())
                     and 'oxidizer_mass_flow_rate' in list(self.parameters.keys())):
            missing_parameters.remove('strain_rate')

        if missing_parameters:
            logger.error('ERROR ! Arguments missing for the case to be launched')
            logger.error('Required arguments are : ' + ', '.join(kwdict.required_args[self.reactor_type]))
            logger.error(', '.join(missing_parameters) + ' is/are missing')
            quit()

        # Unpacking parameters
        self.unpacking()

        # Defining the number of samples
        if 'nsamples' in self.parameters:
            self.nsamples = self.parameters['nsamples']
        else:
            self.nsamples = 50

        # Parameters related to data storing
        self.stored_data = {}

    def copy(self, new_name=None):
        """Copy the class :func:`~ARCANE.mechanisms.Mechanism` object.

        Parameters
        ----------
        new_name :
            name of the new case (Default value = None)

        Returns
        -------
        New class :func:`~ARCANE.mechanisms.Mechanism` object

        """

        mechanism = self.mechanism
        mechanism_copy = self.mechanism.copy()
        if self.mechanism.skeletal:
            skeletal = self.mechanism.skeletal
            skeletal_mechanism_copy = self.mechanism.skeletal.copy()
        self.mechanism = None

        # Re-initialise data
        self.stored_data = {}
        self.reduced_data = {}

        duplicate_case = deepcopy(self)
        if new_name:
            duplicate_case.myid = new_name
            database.load_solutions_names(self, mechanism)
        duplicate_case.mechanism = mechanism_copy
        self.mechanism = mechanism

        if self.mechanism.skeletal:
            duplicate_case.mechanism.skeletal = skeletal_mechanism_copy
            self.mechanism.skeletal = skeletal

        return duplicate_case

    def access_cantera_object(self, mechanism=None):
        """Retrives the Cantera flame object associated to the case.
        
         Parameters
        ----------
        mechanism :
             (Default value = None)

        Returns
        -------
        Cantera flame object

        """

        if not mechanism:
            mechanism = self.mechanism

        ctmech = mechanism.ctmech

        if self.reactor_type in kwdict.reactor_function_name:
            function_name = kwdict.reactor_function_name[self.reactor_type]
            flame_object = getattr(ct, function_name)(ctmech)
        else:
            logger.error('This reactor does not exist. You can only access \
                    premixed flames, counterflow diffusion flames, flamelets \
                    and burners.')
            flame_object = None

        return flame_object

    def unpacking(self):
        """Unpacks the different variables and creates the case id.
        
 
        """

        case_id = self.reactor_type

        for parameter in self.parameters:
            setattr(self, parameter, self.parameters[parameter])

        # Pressure
        if 'pressure' in self.parameters:
            case_id += 'p' + '{:04d}'.format(int(self.pressure / 1000))
        else:
            self.pressure = 1e5

        # Temperature
        if 'temperature' in self.parameters:
            case_id += 'T' + '{:04d}'.format(int(self.temperature))
        #else:
        #    self.temperature = 300

        if 'fuel_temperature' in self.parameters:
            case_id += 'Tfuel' + '{:04d}'.format(int(self.fuel_temperature))

        if 'oxidizer_temperature' in self.parameters:
            case_id += 'Toxi' + '{:04d}'.format(int(self.oxidizer_temperature))

        # Mixture composition
        if 'fuel' in self.parameters:
            # Set up fuel composition for this case
            self.fuel = tools.set_composition_from_string(self.parameters['fuel'],
                                                          self.mechanism.ctmech)
            self.parameters['fuel'] = self.fuel
            fuel_flag = True
        else:
            fuel_flag = False

        if 'oxidizer' in self.parameters:
            # Set up oxidizer composition for this case
            self.oxidizer = tools.set_composition_from_string(self.parameters['oxidizer'],
                                                              self.mechanism.ctmech)
            self.parameters['oxidizer'] = self.oxidizer
            oxidizer_flag = True
        else:
            oxidizer_flag = False

        if 'fuel_air_ratio' in self.parameters and fuel_flag and oxidizer_flag:
            phi = tools.set_phi_from_fuel_air_ratio(self.parameters['fuel_air_ratio'], self.fuel, self.oxidizer,
                                                    self.mechanism.ctmech)
            self.parameters['phi'] = phi

        compo_flag = False
        if 'composition' in self.parameters or 'mass_composition' in self.parameters:
            compo_flag = True
            if 'mass_composition' in self.parameters:
                mass_composition = True
                self.parameters['composition'] = self.parameters['mass_composition']
            else:
                mass_composition = False
            # Set up fuel composition for this case
            self.composition = tools.set_composition_from_string(self.parameters['composition'], self.mechanism.ctmech,
                                                                 dict_in_mass=mass_composition)
            self.parameters['composition'] = self.composition
            if fuel_flag and oxidizer_flag:
                if 'new_phi' in self.parameters and self.parameters['new_phi']:
                    self.phi = tools.set_new_phi_from_composition(self.composition.copy(), self.mechanism.ctmech)
                else:
                    self.phi = tools.set_phi_from_composition(self.composition.copy(), self.fuel.copy(),
                                                              self.oxidizer.copy(), self.mechanism.ctmech)

                self.fuel_air_ratio = tools.set_fuel_air_ratio_from_phi(self.phi, self.fuel.copy(), self.oxidizer.copy(),
                                                                        self.mechanism.ctmech)

        elif 'phi' in self.parameters and fuel_flag and oxidizer_flag:
            compo_flag = True
            if 'new_phi' in self.parameters and self.parameters['new_phi']:
                self.composition = tools.set_composition_from_new_phi(self.parameters['phi'], self.fuel.copy(),
                                                                      self.oxidizer.copy(), self.mechanism.ctmech)
                self.phi = tools.set_new_phi_from_composition(self.composition.copy(), self.mechanism.ctmech)
            else:
                self.composition = tools.set_composition_from_phi(self.parameters['phi'], self.fuel.copy(),
                                                                  self.oxidizer.copy(), self.mechanism.ctmech)

                self.phi = self.parameters['phi']

            self.fuel_air_ratio = tools.set_fuel_air_ratio_from_phi(self.phi, self.fuel.copy(), self.oxidizer.copy(),
                                                                    self.mechanism.ctmech)

        if compo_flag:
            if fuel_flag and oxidizer_flag:
                case_id += 'phi' + '{:03{f}}'.format(self.phi * 100 if self.phi > 1. else int(self.phi*100),
                                                     f='g' if self.phi > 1. else 'd')
            else:
                case_id += 'Compo'
                for spec in self.composition:
                    if int(self.composition[spec] * 100) > 0:
                        case_id += spec + '%{:03d}'.format(int(self.composition[spec] * 100))

        if fuel_flag and oxidizer_flag:
            # Fuel composition
            case_id += 'Fuel'
            for f in self.fuel:
                case_id += f + '%{:03d}'.format(int(self.fuel[f] * 100))

            # Oxidizer composition
            case_id += 'Oxi'
            for o in self.oxidizer:
                case_id += o + '%{:03d}'.format(int(self.oxidizer[o] * 100))

        # Strain rate
        if 'strain_rate' in self.parameters:
            case_id += 'a' + '{:09d}'.format(int(self.strain_rate))

        elif 'fuel_mass_flow_rate' in self.parameters and 'oxidizer_mass_flow_rate' in self.parameters:
            self.fuel_mass_flow_rate = self.parameters['fuel_mass_flow_rate']
            self.oxidier_mass_flow_rate = self.parameters['oxidizer_mass_flow_rate']
            case_id += 'fmdot' + '{:03d}'.format(int(self.fuel_mass_flow_rate * 10))
            case_id += 'omdot' + '{:03d}'.format(int(self.oxidizer_mass_flow_rate * 10))

        elif 'mass_flow_rate' in self.parameters:
            self.mass_flow_rate = self.parameters['mass_flow_rate']
            case_id += 'mdot' + '{:03d}'.format(int(self.mass_flow_rate * 1000))

        elif 'fuel_velocity' in self.parameters and 'oxidizer_velocity' in self.parameters:
            self.fuel_velocity = self.parameters['fuel_velocity']
            self.oxidier_velocity = self.parameters['oxidizer_velocity']
            case_id += 'fv' + '{:03d}'.format(int(self.fuel_velocity * 10))
            case_id += 'ov' + '{:03d}'.format(int(self.oxidizer_velocity * 10))

        elif 'velocity' in self.parameters:
            self.velocity = self.parameters['velocity']
            case_id += 'v' + '{:03d}'.format(int(self.velocity * 1000))

        if 'scalar_dissipation_rate' in self.parameters:
            case_id += 'chist' + '{:09d}'.format(int(self.scalar_dissipation_rate * 1000))

        if 'soot_sections' in self.parameters:
            self.soot_id = 'SOOT' + '{:02d}'.format(self.parameters['soot_sections'])
        if 'soot_precursors' in self.parameters:
            self.soot_id += 'PREC'
            for precursor in self.parameters['soot_precursors']:
                    self.soot_id += precursor
            case_id += self.soot_id

        if 'width' in self.parameters:
            case_id += 'width' + '{:d}'.format(int(self.parameters['width']*100))

        if 'name_addon' in self.parameters:
            case_id += '$' + self.parameters['name_addon']

        # Setting case id
        if 'myid' not in self.parameters:
            self.myid = case_id

        # Setting nickname for plotting
        if 'nickname' not in self.parameters:
            self.nickname = self.myid

        # Targets
        if 'targets' in self.parameters:
            targets = self.targets.copy()
            intersection = list(set(kwdict.names['HR']).intersection(self.targets))
            if intersection:
                targets.remove(intersection[0])

            targets = tools.sort_species_by_solution(targets, self.mechanism.ctmech)

            if intersection:
                targets.append('HeatRelease')

            self.targets = targets

            self.parameters['targets'] = self.targets

        # Error directory
        if 'error_dict' in self.parameters:
            self.check_error_dict()
            self.parameters['error_dict'] = self.error_dict

    def change_case_input_parameter(self, parameter_names, parameter_values):
        """Changing the input parameters of the Case an reloading it.

        Parameters
        ----------
        parameter_names :
            Name of the parameter to change
        parameter_values :
            Value to be changed to
        
        """
        if not type(parameter_names) == list:
            parameter_names = [parameter_names]

        if not type(parameter_values) == list:
            parameter_values = [parameter_values]

        # Changing the parameter
        for name, value in zip(parameter_names, parameter_values):
            self.parameters[name] = value

        # Unpacking parameters
        self.unpacking()

        # Defining the number of samples
        if 'nsamples' in self.parameters:
            self.nsamples = self.parameters['nsamples']
        else:
            self.nsamples = 50

        # Parameters related to data storing
        self.stored_data = {}

        # Instantiating solution name
        database.load_solutions_names(self, mechanism=self.mechanism)

    def check_error_dict(self):
        """
        Verifying the validity of the specified error.
        
        :return: None
        
        """

        possible_keywords = [y for x in kwdict.names.values() for y in x if x not in ['Y', 'X', 'c']]

        for key in self.error_dict:

            if not callable(key):

                quantity, method, type_of_error, location_of_error = error.extract_quantity_method(key)

                if quantity not in possible_keywords:
                    logger.error("WARNING: Error not valid, given keyword: " + key)
                    logger.error('The keywords are :' + str(possible_keywords))
                    quit()

            if not type(self.error_dict[key]) is float and not type(self.error_dict[key]) is int:
                logger.error('ERROR: Non-valid definition of the error, dictionary value must be a float')
                quit()

    def create_profile_from_xml(self, mechanism, solution_name=None, copy_file=True):
        """
        Create a profile file from an existing cantera flame solution
        if a solution_name is given, the file under this path will be restored
        otherwise it looks for the solution with the case id

        mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
        solution_name : str
            file name of the existing cantera solution
        copy_file : bool
            if True, the input file under solution_name is copied as the default xml in the database

        :return:None

        Created: 22/04/08 [TL]
        """

        ctmech = mechanism.ctmech

        if self.reactor_type == 'C1DP':
            flame_object = ct.FreeFlame(ctmech)

        elif self.reactor_type == 'C1DCF':
            flame_object = ct.CounterflowDiffusionFlame(ctmech)

        elif self.reactor_type == 'C1DF':
            flame_object = ct.Flamelet(ctmech)

        elif self.reactor_type == 'C1DB':
            flame_object = ct.BurnerFlame(ctmech)

        else:
            logger.error('This reactor does not exist. You can only rewrite profiles for \
                    premixed flames, counterflow diffusion flames, flamelets \
                    and burners.')

        database.load_solutions_names(self, mechanism)
        if not solution_name:
            solution_name = self.xml_solution
        else:
            subprocess.call(f"cp {solution_name} {self.xml_solution}", shell=True)

        flame_object.restore(solution_name)

        # Save solution
        self.position = flame_object.flame.grid
        self.ndata = len(self.position)
        self.nx = len(self.position)
        self.axis = self.position
        self.T = flame_object.T
        self.P = flame_object.P * np.ones(self.ndata)
        self.Y = flame_object.Y
        self.Y = np.transpose(self.Y)
        self.X = flame_object.X
        if self.reactor_type not in kwdict.reactor_labels['C1DF']:
            self.u = flame_object.u
        self.rho = flame_object.density_mass
        self.nu = flame_object.viscosity
        self.concentrations = flame_object.concentrations
        self.concentrations = np.transpose(self.concentrations)
        self.HR = []
        for index in range(self.nx):
            ctmech.TPX = self.T[index], self.P[index], self.X[:, index]
            self.HR.append(- np.dot(ctmech.net_production_rates, ctmech.partial_molar_enthalpies))
        if hasattr(self, 'soot_sections') and self.soot_sections:
            self.soot_fv = flame_object.soot_fv()
            self.soot_Np = flame_object.soot_Np()


        # Very important part for correctness
        # Storing data
        if len(self.axis) > 0:
            # Corresponding data
            data = [self.axis, self.T, self.P]
            for i in range(len(self.Y[0, :])):
                data.append(self.Y[:, i])
            data += [self.HR]
            if not self.reactor_type in kwdict.reactor_labels['C1DF']:
                data += [self.u]
            data += [self.rho, self.nu]
            if hasattr(self, 'soot_sections') and self.soot_sections:
                data += [self.soot_fv, self.soot_Np]

        # List of data names
        grid_string = kwdict.reactor_grid_name[kwdict.get_reactor_name(self.reactor_type, short=True)]

        data_names = [grid_string, 'Temperature', 'Pressure']
        for spec_name in ctmech.species_names:
            data_names.append(spec_name)
        data_names += ['HeatRelease']

        if not self.reactor_type in kwdict.reactor_labels['C1DF']:
            data_names += ['Velocity']

        data_names += ['Density', 'Viscosity']

        if hasattr(self, 'soot_sections') and self.soot_sections:
            data_names += ['SootVolumeFraction', 'SootNumberDensity']

        # Creating a dictionary
        reduced_data_dict = dict(zip(data_names, data))
        reduced_data_dict['Grid'] = data[0]

        # Dumping data
        database.write_profiles(reduced_data_dict, self, mechanism)

        profile_name = self.solution
        logger.info(f"Profile created as {profile_name}")

        return

    def extract_quantity(self, quantity, mechanism=None):
        """Extracts a scalar value from a profile given the quantity and the method.

        Parameters
        ----------
        quantity :
            Quantity to extract (Laminar flame speed, Temperature max, HR, ...)
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)

        Returns
        -------
        scalar value

        """

        if not mechanism:
            mechanism = self.mechanism

        database.load_solutions_names(self, mechanism)

        if not database.check_solution_existence(self, mechanism):
            return np.nan

        # Extracting quantity and method from string
        quantity, method, type_of_error, location_of_error = error.extract_quantity_method(quantity)

        if location_of_error:
            real_quantity = quantity
            quantity, method, _, _ = error.extract_quantity_method(location_of_error)

        # Take the grid if necessary
        if method in kwdict.methods['int'] + kwdict.methods['dist']:
            x_profile = error.define_error_quantity('grid', self, mechanism)

        # Take the value of the quantity for the reference and the current mechanism
        quantity_profile = error.define_error_quantity(quantity, self, mechanism)

        # Convert them into arrays
        quantity_profile = np.array(quantity_profile)
        quantity_profile = np.nan_to_num(quantity_profile)

        # Minimum of the curve
        if method in kwdict.methods['min']:
            scalar = np.min(quantity_profile)

        # Maximum of the curve
        elif method in kwdict.methods['max']:
            scalar = np.max(quantity_profile)

        # Beginning of the curve
        elif method in kwdict.methods['init']:
            scalar = quantity_profile[0]

        # End of the curve
        elif method in kwdict.methods['end']:
            scalar = quantity_profile[-1]

        # Average of the curve
        elif method in kwdict.methods['mean']:
            scalar = np.mean(quantity_profile)

        # Integral of the curve
        elif method in kwdict.methods['int']:
            scalar = np.trapz(quantity_profile, x_profile)

        else:
            scalar = quantity_profile

        if location_of_error:
            index_of_error = list(quantity_profile).index(scalar)

            real_value = error.define_error_quantity(real_quantity, self, mechanism)

            scalar = real_value[index_of_error]

        return scalar

    def extract_profile_old(self, mechanism):
        """
        Extract data from case simulation results.
        Not to be used, should be removed when path_analysis and hr_analysis are updated to the new sampling method.

        :param mechanism: Mechanism object

        :return data: database

        Created: 11/14/17 [PP]
        Last modified: 11/14/17 [PP]
            - TODO : remove the routine when analysis.py is updated to the new sampling method
        """

        # Instantiating solution name
        database.load_solutions_names(self, mechanism)

        # Go through profile file and extract nsamples states (class attribute)
        if sys.version_info[0] < 3:
            infile = open(self.solution, 'rb')
        else:
            infile = open(self.solution, 'r', newline='', encoding='utf8')

        with infile as csvfile:
            data = np.loadtxt(csvfile, skiprows=2, delimiter=' ')  # , quotechar='|')

        return data

    def names_dictionary_old(self, mechanism):
        """
        Extracts the variables from the file header
        Not to be used, should be removed when path_analysis and hr_analysis are updated to the new sampling method.

        :param mechanism: Mechanism object

        :return: dictionary of name/index correspondence

        Last modified: 27/10/22 [JW]
            - TODO : remove the routine when analysis.py is updated to the new sampling method
        """

        # Instantiating solution name
        database.load_solutions_names(self, mechanism)
        print(self.solution)

        # Go through profile file and extract nsamples states (class attribute)
        if sys.version_info[0] < 3:
            infile = open(self.solution, 'rb')
        else:
            infile = open(self.solution, 'r', newline='', encoding='utf8')

        with infile as csvfile:

            reader = csv.reader(csvfile, delimiter=' ')
            line = next(reader)
            #line = next(reader)
            data_names = line

        # Creating a dictionary
        data_dict = dict(zip(data_names, list(range(len(data_names)))))
        data_dict['Grid'] = 0

        # Duplicates for easier coding
        data_dict['grid'] = 0
        data_dict['T'] = data_dict['Temperature']
        data_dict['P'] = data_dict['Pressure']
        data_dict['HR'] = data_dict['HeatRelease']
        data_dict['HR'] = data_dict['HeatRelease']
        data_names_dict = data_dict.copy()

        return data_names_dict

    def rewrite_solution_from_profile(self, mechanism=None):
        """
        Restores a Flame object from the profile file.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
        flame_object :
            flame that will be computed

        """
        if not mechanism:
            mechanism = self.mechanism
        ctmech = mechanism.ctmech

        # Initialising ctmech to avoid legacy values
        ctmech.TPX = self.temperature, self.pressure, self.composition

        # Loading the data_dict that reads the profile
        data_dict = self.data_dict(mechanism)

        # Retrieving grid
        init_grid = data_dict['Grid']
        width = init_grid[-1]
        flame_grid = np.array(init_grid) / width  # Between 0 and 1 for the set_profile function

        if self.reactor_type == 'C1DP':
            flame_object = ct.FreeFlame(ctmech, grid=init_grid)

        elif self.reactor_type == 'C1DCF':
            flame_object = ct.CounterflowDiffusionFlame(ctmech, grid=init_grid)

        elif self.reactor_type == 'C1DF':
            flame_object = ct.Flamelet(ctmech, grid=init_grid)

        elif self.reactor_type == 'C1DB':
            flame_object = ct.BurnerFlame(ctmech, grid=init_grid)

        else:
            logger.error('This reactor does not exist. You can only rewrite solutions for \
                    premixed flames, counterflow diffusion flames, flamelets \
                    and burners.')

        flame_object.set_profile('T', flame_grid, data_dict['T'])
        flame_object.set_profile('u', flame_grid, data_dict['Velocity'])
        for species in mechanism.species_names:
            flame_object.set_profile(species, flame_grid, data_dict[species])

        loglevel = 1
        flame_object.solve(loglevel=loglevel)
        solution_path, solution_name = database.xml_solution_file_and_name(self.xml_solution)
        flame_object.save(solution_path, name=solution_name, loglevel=loglevel)

        logger.info(f"{self.xml_solution} successfully restored from profile")

    def create_flame_object(self, mechanism):
        """Creates the flame object according to the reactor type.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        `Cantera` Flame object

        """

        ctmech = mechanism.ctmech

        if self.reactor_type == 'C1DP':
            # Creating freely propagating premixed flame
            if hasattr(self, 'width') and hasattr(self, 'initial_grid'):
                if self.width != np.max(self.initial_grid) - np.min(self.initial_grid):
                    logger.error('Error : you have both width and initial_grid set !')
                flame_object = ct.FreeFlame(ctmech, grid=self.initial_grid)
            elif hasattr(self, 'width'):
                flame_object = ct.FreeFlame(ctmech, grid=None, width=self.width)
            elif hasattr(self, 'initial_grid'):
                flame_object = ct.FreeFlame(ctmech, grid=self.initial_grid)
                self.width = np.max(self.initial_grid) - np.min(self.initial_grid)
            else:
                flame_object = ct.FreeFlame(ctmech, width=0.02)
                self.width = 0.02

        elif self.reactor_type == 'C1DCF':
            if hasattr(self, 'width') and hasattr(self, 'initial_grid'):
                if self.width != np.max(self.initial_grid) - np.min(self.initial_grid):
                    logger.error('Error : you have both width and initial_grid set !')
                flame_object = ct.CounterflowDiffusionFlame(ctmech, grid=self.initial_grid)
            elif hasattr(self, 'width'):
                flame_object = ct.CounterflowDiffusionFlame(ctmech, grid=None, width=self.width)
            elif hasattr(self, 'initial_grid'):
                flame_object = ct.CounterflowDiffusionFlame(ctmech, grid=self.initial_grid)
                self.width = np.max(self.initial_grid) - np.min(self.initial_grid)
            else:
                self.initial_grid = np.linspace(0, 0.01, 10)
                flame_object = ct.CounterflowDiffusionFlame(ctmech, grid=self.initial_grid)
                self.width = 0.01

        elif self.reactor_type == 'C1DF':
            if hasattr(self, 'initial_grid'):
                flame_object = ct.Flamelet(ctmech, grid=self.initial_grid)
            else:
                flame_object = ct.Flamelet(ctmech)

        elif self.reactor_type == 'C1DB':
            if hasattr(self, 'width') and hasattr(self, 'initial_grid'):
                if self.width != np.max(self.initial_grid) - np.min(self.initial_grid):
                    logger.error('Error : you have both width and initial_grid set !')
                flame_object = ct.BurnerFlame(ctmech, grid=self.initial_grid)
            elif hasattr(self, 'width'):
                flame_object = ct.BurnerFlame(ctmech, grid=None, width=self.width)
            elif hasattr(self, 'initial_grid'):
                flame_object = ct.BurnerFlame(ctmech, grid=self.initial_grid)
                self.width = np.max(self.initial_grid) - np.min(self.initial_grid)
            else:
                flame_object = ct.BurnerFlame(ctmech, width=0.02)
                self.width = 0.02

        else:
            logger.error('This reactor does not exist. You can only create a flame object for \
                    premixed flames, counterflow diffusion flames, flamelets \
                    and burners.')

        return flame_object

    def clean_solutions(self, mechanism=None):
        """Looks for already computed cases and removes profiles or solution if they are not correct.
        
        
        Parameters
        ----------
        mechanism :
             (Default value = None)

        
        """
        solution_existence = database.check_solution_existence(self, mechanism)

        if os.path.isfile(f"{self.xml_solution}") and not solution_existence:
            subprocess.call(f"rm -rf {self.xml_solution}", shell=True)
            logger.info(f"Solution removed")
        elif not os.path.isfile(f"{self.xml_solution}") and solution_existence:
            self.rewrite_solution_from_profile()
            logger.info(f"XML solution rewritten")
        elif solution_existence and os.path.isfile(f"{self.xml_solution}"):
            data_dict = self.data_dict(mechanism)
            if len(data_dict['Grid']) < 20:
                subprocess.call(f"rm -rf {self.xml_solution} {self.solution}", shell=True)
                logger.info(f"Data removed")
            logger.info("Nothing to be done")
        else:
            logger.info("Nothing to be done")

    def restore_flame(self, mechanism, flame_object):
        """Restores a flame from another.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object
        flame_object :
            flame that will be computed
        from_object :
            tells which object will be use for the restore. It can be either "solution" or "profile"

        Returns
        -------
        boolean telling if the resoration is a success

        """

        if self.restore == 'auto':
            self.restore = 'parent'
            auto_restore = True
        else:
            auto_restore = False

        successful_restore = False

        if self.restore == 'self':
            solution_path, solution_name = database.xml_solution_file_and_name(self.xml_solution)
            if os.path.isfile(solution_path):
                try:
                    flame_object.restore(solution_path, name=solution_name, loglevel=display.cantera_restore_loglevel)
                    logger.info('Flame restored from self')
                    successful_restore = True
                except Exception:
                    logger.warning('Self restoration failed')
            else:
                if auto_restore:
                    self.restore = 'parent'

        if self.restore == 'parent':
            # Automatically restores a solution
            # From a parent mechanism if one exists
            if mechanism.parent:
                solution_path = self.xml_solution.replace(mechanism.name, mechanism.parent.name)

                if os.path.isfile(solution_path):
                    solution_path, solution_name = database.xml_solution_file_and_name(solution_path)
                    flame_object.restore(solution_path, name=solution_name, loglevel=display.cantera_restore_loglevel)
                    logger.info('Flame restored from parent: ' + solution_path)
                    successful_restore = True
                else:
                    if auto_restore:
                        self.restore = 'closest'
            else:
                if auto_restore:
                    self.restore = 'closest'

        # Or from the closest case computed
        if self.restore == 'closest':
            closest_case_id = tools.closest_case_id(self, mechanism,solution=True)
            xml_path = self.xml_solution.replace(self.myid, closest_case_id)
            solution_path, solution_name = database.xml_solution_file_and_name(xml_path)
            if os.path.isfile(solution_path) and not closest_case_id == 'not_valid':
                try:
                    flame_object.restore(solution_path, name=solution_name, loglevel=display.cantera_restore_loglevel)
                    successful_restore = True
                except Exception:
                    logger.warning(f'Automatic restoration of closest case failed ({solution_name})')
                    flame_object = self.create_flame_object(mechanism)
                    pseudo_case = tools.case_from_id(closest_case_id, mechanism)
                    logger.info(f"Regenerating the xml solution file of closest case")
                    pseudo_case.rewrite_solution_from_profile()
                    successful_restore = False
                else:
                    logger.info('Flame restored from closest: ' + xml_path)

        # Restore from user input case name
        solution_path, solution_name = database.xml_solution_file_and_name(self.restore)
        if os.path.isfile(self.restore):
            try:
                flame_object.restore(solution_path, name=solution_name, loglevel=display.cantera_restore_loglevel)
                successful_restore = True
            except Exception:
                logger.warning('Restoration from specified file failed')
                flame_object = self.create_flame_object(mechanism)
                successful_restore = False
            else:
                logger.info('Flame restored from specified: ' + self.restore)
        
        return successful_restore


    def restore_non_sooting(self, mechanism, flame_object):
        """
        Restores a sooting flame from its non-sooting counterpart (computes it if it does not exist)

        :param mechanism: Mechanism object
        :param flame_object: flame that will be computed

        :return: None
        """
        successful_restore = False
        # Recover non-sooting case id
        nonsooting_case_id = self.myid.replace(self.soot_id, '')
        nonsooting_xml_solution = self.xml_solution.replace(self.myid, nonsooting_case_id)
        nonsooting_solution = self.solution.replace(self.myid, nonsooting_case_id)
        # Check if the non-sooting flame exists
        if os.path.isfile(nonsooting_xml_solution):
            try:
                flame_object.restore(nonsooting_xml_solution, loglevel=display.cantera_restore_loglevel)
                successful_restore = True
                logger.info('Non-sooting flame restored from %s' % nonsooting_xml_solution)
            except:
                logger.error('Non-sooting flame %s could not be restored' % nonsooting_xml_solution)
        # Compute the non-sooting flame otherwise
        else:
            sections = self.soot_sections
            case_id = self.myid
            xml_solution = self.xml_solution
            solution = self.solution

            self.soot_sections = 0
            self.myid = nonsooting_case_id
            self.xml_solution = nonsooting_xml_solution
            self.solution = nonsooting_solution
            self.run(mechanism, loglevel=self.loglevel, overwrite=self.overwrite, restore=self.restore)
            self.soot_sections = sections
            self.myid = case_id
            self.xml_solution = xml_solution
            self.solution = solution
            if os.path.isfile(nonsooting_xml_solution):
                flame_object.restore(nonsooting_xml_solution, loglevel=display.cantera_restore_loglevel)
                successful_restore = True
                logger.info('Non-sooting flame computed as %s & restored' % nonsooting_xml_solution)
            else:
                logger.error('Non-sooting flame computed as %s but could not be restored' % nonsooting_xml_solution)  
        
        return successful_restore

    def fill(self, mechanism, data_dict, time=False):
        """Filling a case object with data.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object
        data_dict :
            dictionary of data to be converted
        time :
            if True, the axis will be temporal else spatial (Default value = False)

        """

        self.additional_data = []
        self.additional_data_names = []

        # Parsing input data dictionary
        for key in data_dict:
            if key in kwdict.names['Position'] \
                    + kwdict.names['Time'] \
                    + kwdict.names['grid'] \
                    + kwdict.names['Mixture Fraction']:
                self.axis = data_dict[key]
            elif key in kwdict.names['T']:
                self.T = data_dict[key]
            elif key in kwdict.names['P']:
                self.P = data_dict[key]
            elif key in kwdict.names['HR']:
                self.HR = data_dict[key]
            elif key in mechanism.species_names:
                continue
            else:
                self.additional_data.append(data_dict[key])
                self.additional_data_names.append(key.replace(' ', ''))

        self.ndata = len(self.axis)

        Y_stored = np.zeros([self.ndata, mechanism.ns], 'd')
        for index, species in enumerate(mechanism.species_names):
            if species in data_dict:
                Y_stored[:, index] = data_dict[species]
            else:
                for possible_name in kwdict.names[species]:
                    if possible_name in data_dict:
                        Y_stored[:, index] = data_dict[possible_name]
        self.Y = np.array(Y_stored)

        # Very important part for correctness
        # Storing data
        data = []
        if len(self.axis) > 0:
            # Corresponding data
            data = [self.axis, self.T, self.P]
            for i in range(len(self.Y[0, :])):
                data.append(self.Y[:, i])
            data.append(self.HR)

            data += self.additional_data

        # List of data names
        if time:
            axis = 'Time'
        else:
            axis = 'Position'

        data_names = [axis, 'Temperature', 'Pressure']
        for spec_name in mechanism.ctmech.species_names:
            data_names.append(spec_name)
        data_names.append('HeatRelease')

        data_names += self.additional_data_names

        # Creating a dictionary
        data_dict = dict(zip(data_names, data))

        # Dumping data
        database.write_profiles(data_dict, self, mechanism)

        # Storing data_dict
        self.stored_data[mechanism.name] = data_dict

    def equilibrate(self, mechanism=None, set_as_init=False, variables='HP'):
        """Equilibrates the input mixture of the mechanism.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
        set_as_init :
            if True, the equilibrium mixture becomes the case status (Default value = False)
        variables :
            couple of variables on which equilibrating (Default value = 'HP')


        """
        if not mechanism:
            mechanism = self.mechanism

        ctmech = mechanism.ctmech
        ctmech.TPX = self.temperature, self.pressure, self.composition

        ctmech.equilibrate(variables)

        self.equilibrium = variables
        self.temperature_eq = ctmech.T
        self.pressure_eq = ctmech.P
        self.density_eq = ctmech.density
        self.composition_eq = {}
        self.energy_eq = {}
        for index, spec in enumerate(mechanism.species_names):
            self.composition_eq[spec] = ctmech.X[index]
            self.energy_eq[spec] = ctmech.X[index]

        if set_as_init:
            self.temperature = self.temperature_eq
            self.pressure = self.pressure_eq
            self.composition = self.composition_eq

    def run(self, mechanism=None, loglevel=0, overwrite=False, restore='auto', permissive=True):
        """Running the case:

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
        loglevel :
            level of verbose (Default value = 0)
        overwrite :
            if True, overwrites the case (Default value = False)
        restore :
            method on how a previous solution will be restored (Default value = 'auto')
        permissive :
            if True, failed cases will be computed again (Default value = True)

        """

        self.loglevel = loglevel
        self.overwrite = overwrite
        self.restore = restore

        database.load_solutions_names(self, mechanism=mechanism)

        # Loading default mechanism if needed
        if not mechanism:
            mechanism = self.mechanism

        # Test if this specific case has been run already
        if database.check_solution_existence(self, mechanism) and not overwrite:
            logger.debug(f"Solution {self.solution} already exists: this case has been simulated already")
            self.success = True
            # Get number of data from profile file
            database.read_profiles(self, mechanism)
            if self.ndata < 30 and permissive:
                logger.info(f"Stored solution of case {self.solution} failed and will be redone.")
            else:
                return

        # Compiling the mechanism if needed
        if mechanism.f90 and not mechanism.dynamic_lib_opened:
            mechanism.compile()

        logger.info(f"Running case {self.solution}")

        # Checking if species in inlet composition are in mechanism
        if hasattr(self, 'composition'):
            new_composition = self.composition.copy()
            if hasattr(self, 'variable_inlet'):
                if self.variable_inlet:
                    for spec in self.composition:
                        if spec not in mechanism.ctmech.species_names:
                            new_composition.pop(spec)
                    self.composition = new_composition

        # Calling correct function
        data, data_names = getattr(self, 'compute_' + self.reactor_type)(mechanism)

        # # Resetting the dynamic library
        # if mechanism.f90 and mechanism.dynamic_lib_opened:
        #     mechanism.reset()

        # Creating a dictionary
        reduced_data_dict = dict(zip(data_names, data))
        reduced_data_dict['Grid'] = data[0]

        # Dumping data
        database.write_profiles(reduced_data_dict, self, mechanism)

        # Storing data_dict
        self.reduced_data = reduced_data_dict
        self.reduced_data['mechanism'] = mechanism

    def data_dict(self, mechanism=None, reload=False, extend_data=False, expand_names=False, store=True):
        """Creates a dictionary containing the profiles data :

        Parameters
        ----------
        mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
        reload : bool
            if True, reads again the solution from the file (Default value = False)
        extend_data : bool
            if True, the postproc data will be added to the data_dict (Default value = False)
        expand_names : bool
            if True, alternative names will be added to the data_dict (Default value = False)
        store : bool
            if True, the data_dict corresponding to the mechanism is stored in a per mechanism directory
            (Default value = True)

        Returns
        -------
        dictionary of profiles

        """

        if not mechanism:
            mechanism = self.mechanism

        # Computing everything only if needed
        if mechanism.name not in self.stored_data or reload:

            if self.reduced_data and self.reduced_data['mechanism'].path == mechanism.path:
                reduced_data_dict = self.reduced_data
            else:
                reduced_data_dict = database.read_profiles(self, mechanism)

            data_dict = {}
            data_dict['mechanism'] = mechanism

            self.ndata = len(reduced_data_dict[list(reduced_data_dict.keys())[0]])

            # Checking if the data is consistent
            all_species_in_data = True
            for spec in mechanism.species_names:
                if spec not in reduced_data_dict.keys():
                    all_species_in_data = False

            if not all_species_in_data:
                self.run(mechanism, overwrite=True)
                data_dict = self.data_dict(mechanism, reload=reload)
                return data_dict

            # Associating data to all its base names
            for key in reduced_data_dict:
                data_dict[key] = reduced_data_dict[key]
                if key in mechanism.species_names:
                    data_dict['Y_' + key] = reduced_data_dict[key]
                else:
                    data_dict[kwdict.get_base_name(key)] = reduced_data_dict[key]

            Y = []
            for index in range(self.ndata):
                Y_spec = []
                for spec in mechanism.species_names:
                    Y_spec.append(reduced_data_dict[spec][index])
                Y.append(Y_spec)

            data_dict['Y'] = Y

            # Computing the grid interval
            grid_interval = [0]
            grid_interval += [data_dict['grid'][index] - data_dict['grid'][index - 1] for index in range(self.ndata)][1:]

            for key in kwdict.names['grid_interval']:
                data_dict[key] = grid_interval
            if kwdict.get_reactor_name(self.reactor_type, short=True).startswith('C0'):
                for key in kwdict.names['delta_t']:
                    data_dict[key] = grid_interval
            else:
                for key in kwdict.names['delta_x']:
                    data_dict[key] = grid_interval

            if store:
                self.stored_data[mechanism.name] = data_dict

        else:
            data_dict = self.stored_data[mechanism.name]

        # Updating ndata
        self.ndata = len(data_dict['grid'])

        if extend_data:
            data_dict = postproc.extend_data_dict(data_dict)

        if expand_names:
            data_dict = kwdict.expand_names(data_dict)

        return data_dict

    def compute_C0DX(self, mechanism, reactor_type):
        """Specific routine for 0D Isochoric and Isobaric reactors.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object
        reactor_type :
            constant volume or constant pressure ('CODV', 'CODP', 'V', 'P')

        Returns
        -------
        data matrix, header of data matrix

        """

        # Retrieving Solution object
        ctmech = mechanism.ctmech

        # Initialisation of the mixture
        ctmech.TPX = self.temperature, self.pressure, self.composition

        self.success = True

        # Create the batch reactor (UV is default,  no need to add parameters)
        if reactor_type in ['CODV', 'V']:
            r = ct.IdealGasReactor(ctmech)
        elif reactor_type in ['CODP', 'P']:
            r = ct.IdealGasConstPressureReactor(ctmech)
        else:
            logger.error(f"The reactor type {reactor_type} is not coded")
            logger.error(f"Please choose between constant volume (CODV) and constant pressure (CODP)")
            quit()

        # Now create a reactor network consisting of the single batch reactor
        # Reason: the only way to advance reactors in time is through a network
        sim = ct.ReactorNet([r])

        # Reset ndata counter (will serve as indicator of end of array)
        self.ndata = 0

        # Initialize solution
        if not hasattr(self, 'max_time'):
            self.max_time = 10.0

        # Number of steps saved in solution
        if not hasattr(self, 'nout'):
            self.nout = 10000

        if not hasattr(self, 'time_step'):
            self.time_step = None

        if not hasattr(self, 'time_continuation_factor'):
            self.time_continuation_factor = 2

        # Initialisation for the time stepping
        if self.time_step:
            self.nout = int(np.floor(self.max_time / self.time_step) + 1)

        self.axis = np.zeros(2 * self.nout, 'd')
        self.T = np.zeros(2 * self.nout, 'd')
        self.P = np.zeros(2 * self.nout, 'd')
        self.Y = np.zeros([2 * self.nout, ctmech.n_species], 'd')
        self.concentrations = np.zeros([2 * self.nout, ctmech.n_species], 'd')
        self.HR = np.zeros(2 * self.nout, 'd')
        self.rho = np.zeros(2 * self.nout, 'd')

        # Save initial solution
        nstep = 0
        t = 0.0
        dstore = 1
        condition_fulfilled = False
        self.axis[0] = t
        self.T[0] = r.thermo.T
        self.P[0] = r.thermo.P
        self.Y[0, :ctmech.n_species] = r.thermo.Y
        self.concentrations[0, :ctmech.n_species] = r.thermo.concentrations
        self.HR[0] = 0.0
        self.rho[0] = r.thermo.density

        # Block for the conditional stop of the computation
        if not hasattr(self, 'stop'):
            self.stop = 'ignition'

        # Defines if the computations must continue after the condition is fulfilled
        hard_stop = True
        self.final_time = 10.0

        if self.stop == 'ignition':
            condition = 'self.T[self.ndata] > ' + str(self.T[0] + 500)
            hard_stop = False

        elif 'equilibrium' in self.stop:
            # The conditions specifies a stop when the temperature goes
            # above a percentage of the difference between fresh and burnt gases values
            if '%' in self.stop:
                stop_threshold = float(self.stop.replace('% equilibrium', '')) / 100
            else:
                stop_threshold = 0.5

            self.equilibrate()
            stop_threshold = self.T[0] + (self.temperature_eq - self.T[0]) * stop_threshold

            condition = 'self.T[self.ndata] > ' + str(stop_threshold)
            hard_stop = False

        elif self.stop == 'no':
            condition = "True == False"
        else:
            # Parsing the stop parameter
            stop_elements = self.stop.split(' ')

            stop_quantity = stop_elements[0]
            stop_comparator = stop_elements[1]
            stop_threshold = stop_elements[2]

            if '%' in stop_threshold:
                stop_threshold = float(stop_threshold.replace('%', '')) / 100
                relative_threshold = True
            else:
                relative_threshold = False

            if stop_quantity in kwdict.names['Time']:
                condition = 't ' + stop_comparator + ' ' + str(stop_threshold)
                self.final_time = float(stop_threshold)

            elif stop_quantity in mechanism.species_names:

                species_index = mechanism.species_names.index(stop_quantity)
                if relative_threshold:
                    stop_threshold *= self.Y[0, species_index]

                condition = 'self.Y[self.ndata, species_index] ' + stop_comparator + ' ' + str(stop_threshold)

            elif stop_quantity in kwdict.names['T']:

                if relative_threshold:
                    stop_threshold *= self.T[0]

                condition = 'self.T[self.ndata] ' + stop_comparator + ' ' + str(stop_threshold)

            elif stop_quantity in kwdict.names['P']:

                if relative_threshold:
                    stop_threshold *= self.P[0]

                condition = 'self.P[self.ndata] ' + stop_comparator + ' ' + str(stop_threshold)

            elif self.stop != 'ignition' and 'equilibrium' not in self.stop:
                logger.error("The stop parameter of the computation is not a valid one as it should be a species, "
                             "temperature, pressure, ignition or equilibrium\n "
                             "Falling back to the default 'ignition' value.")
                self.stop = 'ignition'
                condition = 'self.T[self.ndata] > ' + str(self.T[0] + 500)
                hard_stop = False

        # Run only if case succeeded
        if self.success:

            # Iterate in time and save solution
            while t < self.final_time and t < self.max_time:

                if self.time_step:
                    t += self.time_step
                    sim.advance(t)

                else:
                    try:
                        # Advance equations
                        t = sim.step()
                    except Exception:
                        logger.warning('WARNING: Solution may not be converged')
                        break

                    # Update step counter
                    nstep += 1
                    if nstep > 1000000:
                        break

                if eval(condition) and not condition_fulfilled:
                    condition_fulfilled = True
                    if hard_stop:
                        self.final_time = t
                    else:
                        self.final_time = t * self.time_continuation_factor

                # Store if needed
                if nstep % dstore == 0:
                    self.ndata += 1
                    self.axis[self.ndata] = t
                    self.T[self.ndata] = r.thermo.T
                    self.P[self.ndata] = r.thermo.P
                    self.Y[self.ndata, :ctmech.n_species] = r.thermo.Y
                    self.HR[self.ndata] = - np.dot(ctmech.net_production_rates, ctmech.partial_molar_enthalpies)
                    self.rho[self.ndata] = r.thermo.density

                    # Reduce for next time if necessary
                    if self.ndata + 1 == 2 * self.nout:
                        self.ndata = 0
                        for i in range(self.nout):
                            ii = 2 * i + 1
                            self.ndata += 1
                            self.axis[self.ndata] = self.axis[ii]
                            self.T[self.ndata] = self.T[ii]
                            self.P[self.ndata] = self.P[ii]
                            self.Y[self.ndata, :ctmech.n_species] = self.Y[ii, :ctmech.n_species]
                            self.concentrations[self.ndata, :ctmech.n_species] = \
                                self.concentrations[ii, :ctmech.n_species]
                            self.HR[self.ndata] = self.HR[ii]
                            self.rho[self.ndata] = self.rho[self.ndata]

                        # Double dstore
                        dstore *= 2

        # Checking for non-converged solution
        if self.ndata < 10 or self.T[-1] == self.T[0]:
            self.success = False
            logger.warning('WARNING: The simulation failed !')

        # Very important part for correctness
        # Storing data
        data = []
        if len(self.axis) > 0:
            # Corresponding data
            data = [self.axis[:self.ndata + 1], self.T[:self.ndata + 1], self.P[:self.ndata + 1]]
            for i in range(len(self.Y[0, :])):
                data.append(self.Y[:self.ndata + 1, i])
            data += [self.HR[:self.ndata + 1], self.rho[:self.ndata + 1]]

        # List of data names
        data_names = ['Time', 'Temperature', 'Pressure']
        for spec_name in ctmech.species_names:
            data_names.append(spec_name)
        data_names += ['HeatRelease', 'Density']

        return data, data_names

    def compute_C0DV(self, mechanism):
        """Specific routine for 0D Isochoric reactor.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        data matrix, header of data matrix

        """

        return self.compute_C0DX(mechanism, 'V')

    def compute_C0DP(self, mechanism):
        """Specific routine for 0D Isobaric reactor.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        data matrix, header of data matrix

        """

        return self.compute_C0DX(mechanism, 'P')

    def compute_C1DX(self, mechanism, reactor_type):
        """General routine for 1D flames.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object
        reactor_type :
            Flame reactor type:
            - Premixed
            - Counter-flow Diffusion
            - Burner
            - Counter-flow Premixed
            - Impiging jet

        Returns
        -------
        ata matrix, header of data matrix

        """

        # Retrieving Solution object
        ctmech = mechanism.ctmech

        # Defaults parameters for computation
        if not hasattr(self, 'ratio'):
            self.ratio = 2
        if not hasattr(self, 'slope'):
            self.slope = 0.05
        if not hasattr(self, 'curve'):
            self.curve = 0.05
        if not hasattr(self, 'prune'):
            self.prune = 0.01
        if not hasattr(self, 'max_jac_age'):
            self.max_jac_age = 10
        if not hasattr(self, 'timestep'):
            self.timestep = 1e-8
        if not hasattr(self, 'timestep_list'):
            self.timestep_list = [10, 20, 80, 100]
        if not hasattr(self, 'max_time_step_count'):
            self.max_time_step_count = 500
        if not hasattr(self, 'tol_ss'):
            self.tol_ss = None
        if not hasattr(self, 'tol_ts'):
            self.tol_ts = None
        if not hasattr(self, 'refine_param'):
            self.refine_param = 'refine'
        if not hasattr(self, 'thick'):
            self.thick = 1.0
        if not hasattr(self, 'energy'):
            self.energy = 'enabled'
        if not hasattr(self, 'max_grid_points'):
            self.max_grid_points = 1000
        if not hasattr(self, 'soret_enabled'):
            self.soret_enabled = False
        if not hasattr(self, 'radiation_enabled'):
            self.radiation_enabled = False
        if not hasattr(self, 'soot_sections'):
            self.soot_sections = 0
        
        if self.soot_sections:
            if reactor_type in kwdict.reactor_labels['C1DF']:
                self.success = False
                logger.error('Error : Soot computation is not compatible with flames in Z space')
                quit()
            if not hasattr(self, 'soot_precursors'):
                self.success = False
                logger.error('Error : Soot computation needs at least one soot precursor')
                quit()

        # Initialize data
        data = []
        self.success = True

        ########################################################################
        ############             Flame specific calls              #############
        ########################################################################

        if reactor_type in kwdict.reactor_labels['C1DP']:
            flame_function = ct.FreeFlame
        elif reactor_type in kwdict.reactor_labels['C1DCFD']:
            flame_function = ct.CounterflowDiffusionFlame
        elif reactor_type in kwdict.reactor_labels['C1DF']:
            flame_function = ct.Flamelet
        elif reactor_type in kwdict.reactor_labels['C1DB']:
            flame_function = ct.BurnerFlame
        elif reactor_type in kwdict.reactor_labels['C1DCFP']:
            flame_function = ct.CounterflowPremixedFlame
        elif reactor_type in kwdict.reactor_labels['C1DIJ']:
            flame_function = ct.ImpingingJet
        else:
            logger.error(f'Reactor type {reactor_type} is unknown.')
            return

        ########################################################################

        # Creating the flame object with the specified initial grid or width
        if reactor_type in kwdict.reactor_labels['C1DF']:
            width = 1.0  # In Z-space! Must always be equal to 1 in this class!
            self.width = width
            if hasattr(self, 'initial_points'):
                self.initial_grid = np.linspace(0, width, self.initial_points)
            else:
                self.initial_grid = None
            flame_object = flame_function(ctmech, grid=self.initial_grid)
        elif hasattr(self, 'width') and hasattr(self, 'initial_grid'):
            if self.width != np.max(self.initial_grid) - np.min(self.initial_grid):
                logger.error('Error : you have both width and initial_grid set !')
            flame_object = flame_function(ctmech, grid=self.initial_grid, sections=self.soot_sections)
        elif hasattr(self, 'width'):
            flame_object = flame_function(ctmech, grid=None, width=self.width, sections=self.soot_sections)
        elif hasattr(self, 'initial_grid'):
            flame_object = flame_function(ctmech, grid=self.initial_grid, sections=self.soot_sections)
            self.width = np.max(self.initial_grid) - np.min(self.initial_grid)
        else:
            width = kwdict.reactor_default_values[reactor_type]['width']
            flame_object = flame_function(ctmech, width=width, sections=self.soot_sections)
            self.width = width

        if self.thick != 1.0:
            flame_object.flame.thick = self.thick

        if self.soot_sections:
            successful_restore = self.restore_non_sooting(mechanism, flame_object)
        elif self.restore:
            successful_restore = self.restore_flame(mechanism, flame_object)
        else:
            successful_restore = False

        if self.soot_sections:
            flame_object.soot_setup(precursors=self.soot_precursors,
                                    radiation=self.radiation_enabled,
                                    loglevel=self.loglevel)
            self.timestep_list = 'default'
            self.timsestep = 'default'
            self.max_jac_age = 'default'
            self.refine_param = 'disabled'

        # Reset ndata counter
        self.ndata = 0

        # Initialize solution
        self.nout = 1
        self.axis = np.zeros(2 * self.nout, 'd')
        self.T = np.zeros(2 * self.nout, 'd')
        self.P = np.zeros(2 * self.nout, 'd')
        self.Y = np.zeros([2 * self.nout, ctmech.n_species], 'd')
        self.HR = np.zeros(2 * self.nout, 'd')
        self.nx = np.zeros(2 * self.nout, 'd')
        self.rho = np.zeros(2 * self.nout, 'd')
        self.nu = np.zeros(2 * self.nout, 'd')
        self.X = np.zeros([2 * self.nout, ctmech.n_species], 'd')
        self.u = np.zeros(2 * self.nout, 'd')
        self.position = np.zeros(2 * self.nout, 'd')
        self.concentrations = np.zeros([2 * self.nout, ctmech.n_species], 'd')


        # Set computation parameters
        if hasattr(self, 'imposed_profile'):
            flame_object.energy_enabled = False
            try:
                flame_object.flame.set_fixed_temp_profile(
                        (self.imposed_profile[:, 0] - self.imposed_profile[0, 0])
                        / (self.imposed_profile[-1, 0] - self.imposed_profile[0, 0])
                        , self.imposed_profile[:, 1])
            except Exception:
                self.success = False
                logger.warning('WARNING: Imposed temperature is badly/not formatted !')
        else:
            flame_object.energy_enabled = True

        flame_object.set_refine_criteria(ratio=self.ratio,
                                         slope=self.slope,
                                         curve=self.curve,
                                         prune=self.prune)

        if not self.max_jac_age == 'default':
            flame_object.set_max_jac_age(self.max_jac_age, self.max_jac_age)

        if not self.timestep == 'default' and not self.timestep_list == 'default':
            flame_object.set_time_step(self.timestep, self.timestep_list)

        if not self.max_time_step_count == 'default':
            flame_object.max_time_step_count = self.max_time_step_count

        if hasattr(self, 'grid_min'):
            flame_object.set_grid_min(self.grid_min)

        flame_object.flame.set_steady_tolerances(default=self.tol_ss)
        flame_object.flame.set_transient_tolerances(default=self.tol_ts)

        flame_object.set_max_grid_points(flame_object.flame, self.max_grid_points)

        ########################################################################
        ############             Flame specific calls              #############
        ########################################################################

        # Setting boundary condition
        if reactor_type in kwdict.reactor_labels['C1DP']:

            flame_object.inlet.X = self.composition
            flame_object.inlet.T = self.temperature
            flame_object.P = self.pressure

        elif reactor_type in kwdict.reactor_labels['C1DCFD']:

            flame_object.fuel_inlet.X = self.fuel
            flame_object.fuel_inlet.T = self.fuel_temperature

            flame_object.oxidizer_inlet.X = self.oxidizer
            flame_object.oxidizer_inlet.T = self.oxidizer_temperature

            flame_object.P = self.pressure

            if hasattr(self, 'fuel_velocity') and hasattr(self, 'oxidizer_velocity'):
                ctmech.TPX = self.fuel_temperature, ctmech.P, self.fuel
                self.fuel_mass_flow_rate = self.fuel_velocity * ctmech.density

                ctmech.TPX = self.oxidizer_temperature, ctmech.P, self.oxidizer
                self.oxidizer_mass_flow_rate = self.oxidizer_velocity * ctmech.density

            elif not hasattr(self, 'fuel_mass_flow_rate') and not hasattr(self, 'oxidizer_mass_flow_rate'):
                if not hasattr(self, 'phi'):
                    self.phi = 1
                self.fuel_mass_flow_rate, self.oxidizer_mass_flow_rate = tools.set_mass_flow_rates(ctmech,
                                                                                                   self.phi,
                                                                                                   self.strain_rate,
                                                                                                   self.fuel,
                                                                                                   self.oxidizer,
                                                                                                   self.fuel_temperature,
                                                                                                   self.oxidizer_temperature,
                                                                                                   self.width)

            flame_object.fuel_inlet.mdot = self.fuel_mass_flow_rate
            flame_object.oxidizer_inlet.mdot = self.oxidizer_mass_flow_rate

        elif reactor_type in kwdict.reactor_labels['C1DF']:

            flame_object.fuel_inlet.X = self.fuel
            flame_object.fuel_inlet.T = self.fuel_temperature

            flame_object.oxidizer_inlet.X = self.oxidizer
            flame_object.oxidizer_inlet.T = self.oxidizer_temperature

            flame_object.ChiSt = self.scalar_dissipation_rate

            flame_object.P = self.pressure

        elif reactor_type in kwdict.reactor_labels['C1DB']:

            flame_object.burner.X = self.composition
            flame_object.burner.T = self.temperature
            flame_object.P = self.pressure

            if hasattr(self, 'velocity'):
                flame_object.burner.mdot = self.velocity * ctmech.density
            elif hasattr(self, 'mass_flow_rate'):
                flame_object.burner.mdot = self.mass_flow_rate
            else:
                self.success = False
                logger.warning('WARNING: Specify either mass flow rate or velocity for burner flames')

        elif reactor_type in kwdict.reactor_labels['C1DCFP']:

            ctmech.TPX = self.temperature, self.pressure, self.composition

            density_fresh = ctmech.density

            self.equilibrate()

            density_burnt = self.density_eq

            fact_rho = (1.0 + density_fresh / density_burnt)

            vel_fresh = self.strain_rate * self.width / fact_rho

            mdot_fresh = vel_fresh * density_fresh
            mdot_burnt = mdot_fresh

            flame_object.reactants.mdot = mdot_fresh
            flame_object.products.mdot = mdot_burnt

        elif reactor_type in kwdict.reactor_labels['C1DIJ']:

            flame_object.inlet.X = self.composition
            flame_object.inlet.T = self.temperature

            flame_object.P = self.pressure

            if hasattr(self, 'velocity'):
                flame_object.inlet.mdot = self.velocity * ctmech.density
            elif hasattr(self, 'mass_flow_rate'):
                flame_object.inlet.mdot = self.mass_flow_rate
            else:
                self.success = False
                logger.warning('WARNING : Specify either mass flow rate or velocity for iminging jet flames')
            
            if hasattr(self, 'surface_temperature'):
                flame_object.surface.T = self.surface_temperature
            else:
                self.sucess = False
                logger.warning('WARNING : Specify impinged surface temperature for impinging jet flames')

        else:
            logger.error(f'Reactor type {reactor_type} is unknown.')
            return

        ########################################################################

        if mechanism.transport == 'Multi':
            flame_object.flame.soret_enabled = self.soret_enabled

        if not successful_restore:
            flame_object.set_initial_guess()

        if self.energy == 'disabled' or self.energy == 'disabled_enabled':
            flame_object.energy_enabled = False
            # Solving flame
            logger.debug("Energy equation disabled.")
            try:
                flame_object.solve(loglevel=self.loglevel,
                                   refine_grid=self.refine_param)  # only if you have the right Cantera version
            except Exception:
                self.success = False
                logger.warning('WARNING: The simulation without energy equation failed !')

            if self.success:
                solution_path, solution_name = database.xml_solution_file_and_name(self.xml_solution)
                flame_object.save(solution_path, name=solution_name, loglevel=self.loglevel)

        if self.energy == 'enabled' or self.energy == 'disabled_enabled':
            logger.debug("Energy equation enabled.")

            try:
                flame_object.energy_enabled = True
                flame_object.solve(loglevel=self.loglevel,
                                   refine_grid=self.refine_param)  # only if you have the right Cantera version

                solution_path, solution_name = database.xml_solution_file_and_name(self.xml_solution)
                flame_object.save(solution_path, name=solution_name, loglevel=self.loglevel)

            except Exception:
                self.success = False
                logger.warning('WARNING: The simulation with energy equation failed !')

        # Save solution
        self.position = flame_object.flame.grid
        self.ndata = len(self.position)
        self.nx = len(self.position)
        self.axis = self.position
        self.T = flame_object.T
        self.P = flame_object.P * np.ones(self.ndata)
        self.Y = flame_object.Y
        self.Y = np.transpose(self.Y)
        self.X = flame_object.X
        if reactor_type not in kwdict.reactor_labels['C1DF']:
            self.u = flame_object.u
        self.rho = flame_object.density_mass
        self.nu = flame_object.viscosity
        self.concentrations = flame_object.concentrations
        self.concentrations = np.transpose(self.concentrations)
        self.HR = []
        for index in range(self.nx):
            ctmech.TPX = self.T[index], self.P[index], self.X[:, index]
            self.HR.append(- np.dot(ctmech.net_production_rates, ctmech.partial_molar_enthalpies))
        if self.soot_sections:
            self.soot_fv = flame_object.soot_fv()
            self.soot_Np = flame_object.soot_Np()


        # Very important part for correctness
        # Storing data
        if len(self.axis) > 0:
            # Corresponding data
            data = [self.axis, self.T, self.P]
            for i in range(len(self.Y[0, :])):
                data.append(self.Y[:, i])
            data += [self.HR]
            if not reactor_type in kwdict.reactor_labels['C1DF']:
                data += [self.u]
            data += [self.rho, self.nu]
            if self.soot_sections:
                data += [self.soot_fv, self.soot_Np]

        # List of data names
        grid_string = kwdict.reactor_grid_name[kwdict.get_reactor_name(reactor_type, short=True)]

        data_names = [grid_string, 'Temperature', 'Pressure']
        for spec_name in ctmech.species_names:
            data_names.append(spec_name)
        data_names += ['HeatRelease']

        if not reactor_type in kwdict.reactor_labels['C1DF']:
            data_names += ['Velocity']

        data_names += ['Density', 'Viscosity']

        if self.soot_sections:
            data_names += ['SootVolumeFraction', 'SootNumberDensity']

        return data, data_names

    def compute_C1DP(self, mechanism):
        """Specific routine for 1D Premixed laminar flame.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        data matrix, header of data matrix

        """

        return self.compute_C1DX(mechanism, 'C1DP')

    def compute_C1DCFD(self, mechanism):
        """
        Specific routine for 1D CounterFlow Diffusion Flame.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        data matrix, header of data matrix

        """

        return self.compute_C1DX(mechanism, 'C1DCFD')

    def compute_C1DF(self, mechanism):
        """Specific routine for 1D CounterFlow Diffusion Flame in Z space.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        data matrix, header of data matrix

        """

        return self.compute_C1DX(mechanism, 'C1DF')

    def compute_C1DB(self, mechanism):
        """Specific routine for 1D Burner laminar flame.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        data matrix, header of data matrix

        """

        return self.compute_C1DX(mechanism, 'C1DB')

    def compute_C1DCFP(self, mechanism):
        """
        Specific routine for 1D Burner laminar flame.

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        data matrix, header of data matrix

        """

        return self.compute_C1DX(mechanism, 'C1DCFP')

    def compute_C1DIJ(self, mechanism):
        """
        Specific routine for 1D Impinging Jet laminar flame

        Parameters
        ----------
        mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------
        data matrix, header of data matrix

        """

        return self.compute_C1DX(mechanism, 'C1DIJ')


    def compute_CXC(self, mechanism):
        """
        Specific routine for custom functions

        Parameters
        ----------

        mechanism: class :func:`~ARCANE.mechanisms.Mechanism` object
            class :func:`~ARCANE.mechanisms.Mechanism` object

        Returns
        -------

        data matrix, header of data matrix

        """

        self.ndata = -1
        self.success = True

        function_arguments = self.parameters.copy()

        # Discards parameters that are not related to the function itself
        function_arguments.pop('function')
        if 'targets' in function_arguments:
            function_arguments.pop('targets')
        if 'error_dict' in function_arguments:
            function_arguments.pop('error_dict')

        data, data_names = self.function(mechanism, **function_arguments)
        self.ndata = len(data[0])

        return data, data_names
