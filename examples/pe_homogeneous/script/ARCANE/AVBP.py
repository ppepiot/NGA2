"""Routines for direct integration into AVBP

Author: QC (2018/12/19)
"""

import ARCANE.cases as cases
import ARCANE.mechanisms as mechanisms
import ARCANE.error as error
import ARCANE.database as database
import ARCANE.graphs as graphs
import ARCANE.tools as tools
import ARCANE.display as display
import ARCANE.chem_data as chem_data

import cantera as ct

import scipy.optimize as opt
import numpy as np

import subprocess
import os

from string import Template

logger = display.Logger()
logger.set_log('logCompute')

graphs = graphs.Graph()


def dump_AVBP_inputs(cases_list, mechanism, root=None, output_dir="AVBP_input", overwrite=False,
                     header='Analytically Reduced Chemistry', mixture_name=None,
                     semi_implicit=False, rrate=True, exponential=False,
                     save_init=True, method='cases', author='Unknown Author', since='', note='',
                     fuel_name='FUEL', liquid=[], liquid_properties={}, fuel_is_liquid=False, restore_cases='auto'):

    """
    Writes species_database.dat and mixture_database.dat files

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object
    root : class :func:`~ARCANE.mechanisms.Mechanism` object
        transport comparison will be done with the specified Mechanism object (Default value = None)
    output_dir : str
        directory where the files will be written (Default value = "AVBP_input")
    overwrite : bool
        if False, overwrites any existing file else will append new data (Default value = False)
    header : str
        mixture_database.dat header (also used as details for the f90 header) (Default value = 'Analytically Reduced Chemistry')
    mixture_name : str
        name of the mixture (Default value = None)
    semi_implicit : bool
        list of species that will be semi-implicited (all species if True) (Default value = False), following the method described in T. Jaravel Phd Thesis, available at <https://cerfacs.fr/wp-content/uploads/2016/06/THESE_CFD_JARAVEL.pdf>
    rrate : bool
        if True, writes rrate routine for AVBP (Default value = True)
    exponential : bool
        if True uses the exponential formulation for production rates (Default value = False)
    save_init : bool
        if True, saves initial composition for direct copy/paste in solutbound file (Default value = True)
    method : str
        method for optimization ('basic', 'cases' or 'optimization') (Default value = 'cases')
    author : str
        name of the mechanism author (for AVBP name of the mixture) (Default value = 'Unknown Author')
    since : str
        current version of AVBP (useful to fill the f90 header) (Default value = '')
    note : str
        note on reduction (useful to fill the f90 header) (Default value = '')
    fuel_name : str
        name of the fuel (for AVBP name of the mixture) (Default value = 'FUEL')
    liquid : list
        name of the fuels that have to be set as liquid (Default value = [])
    liquid_properties : dict
        dictionary with liquid properties in the form
        {'species': {'liquid_density': 0,
        'boiling_temperature': 0
        'critical_temperature': 0
        'critical_pressure': 0
        'critical_volume': 0
        'viscosity_params': [0, 0, 0, 0]}} (Default value = {})
    fuel_is_liquid: bool
        if True, all the fuel species are set as liquid
    restore_cases: bool
        defines the method for the cases restoration


    Created: 18/12/15 [QC]

    """

    # Create output directory if no already there (do not overwrite)
    database.create_dir(output_dir, overwrite=False)

    if type(cases_list) != list:
        cases_list = [cases_list]

    # Extracting fuel species
    fuel_species = list(cases_list[0].fuel.keys())
    if fuel_name == 'FUEL':
        fuel_name = fuel_species[0]
        warning_text = '\nWARNING : You should specify '
        warning_text += "a fuel mixture name with the keyword fuel_name='FUEL' otherwise first fuel is used as fuel name"
        logger.warning(warning_text)

    if not mixture_name:
        warning_text = '\nWARNING : You should specify '
        plural_text = 0
        if author == 'Unknown Author':
            if 'ARCANE_AUTHOR' in os.environ:
                author = os.environ.get('ARCANE_AUTHOR')
            else :
                if plural_text == 1:
                    warning_text += '\nand '
                warning_text += "an author name with the keyword author='Wade Wilson'"
                plural_text += 1

        if plural_text > 0:
            logger.warning(warning_text)

        # Creating initials from full_name
        initials_list = author.split(' ')
        if not len(initials_list) == 1:
            initials = ''
            for part in initials_list:
                initials += part[0]
        else:
            initials = author

        mixture_name = fuel_name + '_' + str(mechanism.ns) \
                       + '_' + str(mechanism.nr + mechanism.nr_reverse) \
                       + '_' + str(mechanism.nqss) + '_' + initials

    # Filtering only 1D cases
    new_cases_list = []
    for case in cases_list:
        if not case.reactor_type.startswith('C0'):
            new_cases_list.append(case)

    if cases_list != new_cases_list:
        logger.info('0D cases have been discarded from optimization as they do not involve transport')

    cases_list = new_cases_list

    # Running the reference cases
    cases.run_cases(cases_list, mechanism)

    # Sets the reference mechanism for comparison
    if not root:
        root = mechanism

    # Selects the optimization method
    if method == 'basic':
        best_combination, state = optimize_on_cases_basic(cases_list, mechanism, root=root, overwrite=overwrite)
    elif method == 'cases':
        best_combination, state = optimize_on_cases(cases_list, mechanism, root=root, overwrite=overwrite)
    elif method == 'none':
        best_combination, state = extract_transport_parameters(cases_list, mechanism)
    else:
        best_combination, state = optimize_transport(cases_list, mechanism, root=root)

    if method != 'none':
        # Computing again the best case
        # Writing corresponding mixture_database.dat
        logger.info('\n-----------')
        logger.info(' Best case ')
        logger.info('-----------')
        transport_error(cases_list, mechanism, best_combination, state,
                        root=root, output_dir=output_dir, save_init=save_init)

    logger.info('\nWriting final mixture_database.dat, species_database.dat and fortran file')

    if fuel_is_liquid:
        for key in cases_list[0].fuel:
            if not key in liquid:
                liquid.append(key)

    mechanism.ctmech.X = cases_list[0].fuel

    fuel_Y = {}
    for key in cases_list[0].fuel:
        fuel_Y[key] = mechanism.ctmech.Y[mechanism.ctmech.species_index(key)]

    mechanism.ctmech.X = cases_list[0].oxidizer
    oxidizer_Y = {}
    for key in cases_list[0].oxidizer:
        oxidizer_Y[key] = mechanism.ctmech.Y[mechanism.ctmech.species_index(key)]

    # Writes corresponding mixture database
    write_mixture_database(mechanism, transport_dict=best_combination, state=state,
                           header=header, mixture_name=mixture_name, parsing=True,
                           fuel=fuel_Y, oxidizer=oxidizer_Y)

    # Writes species database
    write_species_database(mechanism, fuel=cases_list[0].fuel.keys(), liquid=liquid,
                           liquid_properties=liquid_properties)

    # Header filling information
    supplementary_info = {}
    supplementary_info['details'] = header
    supplementary_info['authors'] = author
    supplementary_info['since'] = since
    supplementary_info['note'] = note

    if save_init:
        for case in cases_list:
            AVBP_init(case, mechanism, sol_dir=output_dir + "/sol_can2av")
            write_inlet_for_solutbound(case, mechanism, sol_dir=output_dir + "/solutbound_inlet")
            write_inlet_for_makesolution(case, mechanism, sol_dir=output_dir + "/sol_makesolution")
            write_outlet_for_gas_out(case, mechanism, sol_dir=output_dir + "/sol_gas_out")


    # Writes f90 file
    mechanism.write_f90(file_name=mixture_name, use='AVBP', routine_name=mixture_name,
                        mixture_name=mixture_name, fuel=fuel_name, author=author,
                        semi_implicit=semi_implicit, rrate=rrate,
                        exponential=exponential, supplement_dict=supplementary_info)

    if mechanism.f90:
        subprocess.call('mv ' + mixture_name + '.f90 ' + output_dir + '/', shell=True)

    subprocess.call('mv mixture_database.dat ' + output_dir + '/', shell=True)
    subprocess.call('mv species_database.dat ' + output_dir + '/', shell=True)


def transport_extrema(cases_list, mechanism):
    """Computes the maximum and minimum values of the Prandtl and Schmidt numbers on a list of cases

    Parameters
    ----------
    cases_list :
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    transport_dic : dict
        dictionary of extremum values

    """

    min_Sc = [10] * mechanism.ns
    max_Sc = [0] * mechanism.ns
    min_Pr = 10
    max_Pr = 0

    transport_dict = {}

    for index, spec in enumerate(mechanism.ctmech.species_names):
        # Loop on all cases
        for case in cases_list:
            # Extracting the Schmidt profile for one species
            Sc_profile = tools.Schmidt_profile(case, mechanism, spec)
            Sc_profile = [val for val in Sc_profile if not np.isnan(val)]
            min_value = min(Sc_profile)
            max_value = max(Sc_profile)

            # Storing extremum values
            if min_value < min_Sc[index]:
                min_Sc[index] = min_value

            if max_value > max_Sc[index]:
                max_Sc[index] = max_value

            # Extracting Prandtl profile (only on first iteration)
            if index == 0:
                Pr_profile = tools.Prandtl_profile(case, mechanism)
                Pr_profile = [val for val in Pr_profile if not np.isnan(val)]
                min_value = min(Pr_profile)
                max_value = max(Pr_profile)

                # Storing extremum values
                if min_value < min_Pr:
                    min_Pr = min_value

                if max_value > max_Pr:
                    max_Pr = max_value

        # Storing data in dictionary
        transport_dict['min_Sc_' + spec] = min_Sc[index]
        transport_dict['max_Sc_' + spec] = max_Sc[index]
    transport_dict['min_Pr'] = min_Pr
    transport_dict['max_Pr'] = max_Pr

    return transport_dict


def extract_transport_parameters(cases_list, mechanism):
    """Extracts the Prandtl and Schmidt number from the first case of the list

    Parameters
    ----------
    cases_list :
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    best_combination : dict
        dictionary of the transport values for the first case
    state : list
        list of temperature, pressure and mass_fractions corresponding to the first case


    """

    case = cases_list[0]

    best_combination = {}

    species_names = mechanism.species_names

    transport_variables = species_names + ['Pr']

    # Initializing the transport values with the burnt gases of first case
    if case.reactor_type != 'C1DP':
        index_gas_state = [-1]
    else:
        index_gas_state = [1]

    temperature = case.extract_quantity('T', mechanism=mechanism)[index_gas_state]

    pressure = case.extract_quantity('P', mechanism=mechanism)[index_gas_state]

    mass_fractions = case.extract_quantity('Y', mechanism=mechanism)[index_gas_state]

    mechanism.ctmech.TPY = temperature, pressure, mass_fractions

    schmidt = tools.compute_Schmidt_numbers(mechanism.ctmech)
    prandtl = tools.compute_Prandtl_number(mechanism.ctmech)

    for index, var in enumerate(transport_variables):
        if var != 'Pr':
            best_combination[var] = schmidt[index]
        else:
            best_combination[var] = prandtl

    state = [temperature, pressure, mass_fractions]

    return best_combination, state


def optimize_transport(cases_list, mechanism, root=None, n_steps=10, number_of_iterations=2):
    """Optimizes the Prandtl and Schmidt number on a given list of cases
    by varying individually the values of Prandtl and Schmidt numbers

    Parameters
    ----------
    cases_list :
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    root :
        class :func:`~ARCANE.mechanisms.Mechanism object to be compared to (Default value = None)
    n_steps :
        number of steps in extrema range of the transport number (Default value = 10)
    number_of_iterations :
        number of iteration of the optimization routine (Default value = 2)

    Returns
    -------
    best_combination : dict
        dictionary of optimized transport values
    state : list
        list of temperature, pressure and mass_fractions corresponding to optimized transport values

    """

    if not root:
        root = mechanism

    species_names = mechanism.ctmech.species_names

    # Generating list of n values between the extrema values of transport numbers
    transport_dict = transport_extrema(cases_list, mechanism)

    transport_lists = {}

    for spec in species_names:
        step = (transport_dict['max_Sc_' + spec] - transport_dict['min_Sc_' + spec]) / n_steps
        transport_lists[spec] = np.arange(transport_dict['min_Sc_' + spec],
                                          transport_dict['max_Sc_' + spec] + step,
                                          step)

    step = (transport_dict['max_Pr'] - transport_dict['min_Pr']) / n_steps
    transport_lists['Pr'] = np.arange(transport_dict['min_Pr'],
                                      transport_dict['max_Pr'] + step / 2,
                                      step)

    best_combination, state = optimize_on_cases_basic(cases_list, mechanism, root=root)

    # Storing current state to be used for viscosity
    initial_state = [mechanism.ctmech.T, mechanism.ctmech.P, mechanism.ctmech.Y]

    # Writing dummy mixture_database
    write_mixture_database(mechanism, transport_dict=best_combination)

    # Creating Mechanism object with AVBP transport
    mechanism_AVBP = mechanisms.Mechanism(mechanism.path, transport='AVBP')
    mechanism_AVBP.name += '_AVBP'

    # Sorting species from most to least abundant
    mass_fractions_integrals = {}
    for spec in species_names:
        mass_fractions_integrals[spec] = 1e-60

    for spec in species_names:
        for case in cases_list:
            integral = case.extract_quantity(spec + ' int', mechanism=mechanism)
            if integral > mass_fractions_integrals[spec]:
                mass_fractions_integrals[spec] = integral

    sorted_tuple = sorted(zip(list(mass_fractions_integrals.values()), mass_fractions_integrals.keys()),
                          key=lambda x: x[0])

    species_names = [x[1] for x in sorted_tuple]
    species_names.reverse()

    # Inverting for better results
    transport_variables = ['Pr'] + species_names
    mass_fractions_integrals['Pr'] = 1

    for iteration in range(number_of_iterations):
        if iteration > 0:
            logger.info('\nEntering iteration ' + str(iteration + 1) + '\n')
        for var in transport_variables:
            if var != 'Pr':
                logger.info('\nComputing variations of ' + var + ' Schmidt number\n')
            else:
                logger.info('\nComputing variations of Prandtl number\n')

            # Previous logger level
            previous_level_index = logger.level_index
            logger.set_log(i=previous_level_index + 1, log_name='logCompute')
            logger.set_log(i=previous_level_index + 2, log_name='logReduction')

            # If the value does not change by at least 1%, it want be varied
            variation = (transport_lists[var][-1] - transport_lists[var][0]) / transport_lists[var][0]
            if variation > 0.01 and mass_fractions_integrals[var] > 1e-5:

                # Storing best value
                var_error_square = []
                for step in range(n_steps):

                    # Incrementing one transport parameter
                    best_combination[var] = transport_lists[var][step]

                    # Writing corresponding mixture_database.dat
                    write_mixture_database(mechanism, transport_dict=best_combination)
                    # Re-creating Mechanism object with AVBP transport to update values
                    mechanism_AVBP = mechanisms.Mechanism(mechanism.path, transport='AVBP')
                    mechanism_AVBP.name += '_AVBP'

                    sum_square = 0
                    for case in cases_list:
                        # Loading reference case xml solution
                        database.load_solutions_names(case, mechanism)
                        file_to_restore = case.xml_solution

                        case.run(mechanism_AVBP, restore=file_to_restore, overwrite=True)

                        error_list = error.case_error(case, root, mechanism=mechanism_AVBP)

                        # Computing sum of the squared errors weighed by their tolerances
                        weighed_error_list = [(error / tolerance) ** 2
                                              for error, tolerance in zip(error_list, list(case.error_dict.values()))]
                        sum_square += sum(weighed_error_list)

                    # Storing the errors
                    var_error_square.append(sum_square)

                # Choosing the value with the minimal error for the following
                index_min = var_error_square.index(min(var_error_square))
                best_combination[var] = transport_lists[var][index_min]

            logger.set_log(i=previous_level_index, log_name='logCompute')
            logger.set_log(i=previous_level_index, log_name='logReduction')

    state = [initial_state[0], initial_state[1], initial_state[2]]
    mechanism.ctmech.TPY = initial_state[0], initial_state[1], initial_state[2]

    return best_combination, state


def optimize_on_cases_basic(cases_list, mechanism, root=None, overwrite=False, log=3):
    """Optimization of the Prandtl and Schmidt number done by assessing the error each case produce on itself

    Parameters
    ----------
    cases_list :
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    root :
        class :func:`~ARCANE.mechanisms.Mechanism` object to be compared to (Default value = None)
    overwrite :
        if True, overwrites already computed cases (Default value = False)
    log :
        level of display (Default value = 3)

    Returns
    -------
    best_combination : dict
        dictionary of optimized transport values
    state : list
        list of temperature, pressure and mass_fractions corresponding to optimized transport values

    """

    if not root:
        root = mechanism

    species_names = mechanism.ctmech.species_names

    transport_variables = species_names + ['Pr']

    combination = {}

    var_error_square = []
    store_combination = []

    # Storing states
    store_temperatures = []
    store_pressures = []
    store_mass_fractions = []

    sum_square = 0
    for case_to_run in cases_list:

        # Initializing the transport values with the burnt gases of first case
        if case_to_run.reactor_type != 'C1DP':
            temperatures = list(case_to_run.extract_quantity('T', mechanism=mechanism))
            temperature = max(temperatures)
            index_gas_state = temperatures.index(temperature)
            if type(index_gas_state) == list():
                index_gas_state = index_gas_state[0]
            indices_gas_state = [index_gas_state]
            list_of_states = ['burnt']
        else:
            indices_gas_state = [1, -1]
            list_of_states = ['fresh', 'burnt']

        for index_in_list, index_gas_state in enumerate(indices_gas_state):

            temperature = case_to_run.extract_quantity('T', mechanism=mechanism)[index_gas_state]
            store_temperatures.append(temperature)

            pressure = case_to_run.extract_quantity('P', mechanism=mechanism)[index_gas_state]
            store_pressures.append(pressure)

            mass_fractions = case_to_run.extract_quantity('Y', mechanism=mechanism)[index_gas_state]
            store_mass_fractions.append(mass_fractions)

            mechanism.ctmech.TPY = temperature, pressure, mass_fractions

            schmidt = tools.compute_Schmidt_numbers(mechanism.ctmech)
            prandtl = tools.compute_Prandtl_number(mechanism.ctmech)

            for index, var in enumerate(transport_variables):
                if var != 'Pr':
                    combination[var] = schmidt[index]
                else:
                    combination[var] = prandtl

            # Writing dummy mixture_database
            write_mixture_database(mechanism, transport_dict=combination)

            # Creating Mechanism object with AVBP transport
            mechanism_AVBP = mechanisms.Mechanism(mechanism.path, name=f'{mechanism.name}_AVBP_{case_to_run.myid}'
                                                     f'_{list_of_states[index_in_list]}', transport='AVBP')

            # Loading reference case xml solution
            database.load_solutions_names(case_to_run, mechanism)
            file_to_restore = case_to_run.xml_solution

            case_to_run.run(mechanism_AVBP, restore=file_to_restore, overwrite=overwrite)

            error_list = error.case_error(case_to_run, root, mechanism=mechanism_AVBP)

            # Computing sum of the squared errors weighed by their tolerances
            weighed_error_list = [(error / tolerance) ** 2
                                  for error, tolerance in zip(error_list, list(case_to_run.error_dict.values()))]
            sum_square += sum(weighed_error_list)

            # Storing the errors
            var_error_square.append(sum_square)
            store_combination.append(combination.copy())

    index_best = var_error_square.index(min(var_error_square))

    # Setting the mechanism to the correct state
    state = [store_temperatures[index_best], store_pressures[index_best], store_mass_fractions[index_best]]
    mechanism.ctmech.TPY = store_temperatures[index_best], store_pressures[index_best], store_mass_fractions[index_best]

    best_combination = store_combination[index_best]

    return best_combination, state


def optimize_on_cases(cases_list, mechanism, root=None, overwrite=True):
    """Optimization of the Prandtl and Schmidt number done by assessing the error each case produce on the others

    Parameters
    ----------
    cases_list :
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    root :
        class :func:`~ARCANE.mechanisms.Mechanism` object to be compared to (Default value = None)
    overwrite :
        if True, overwrites already computed cases (Default value = True)

    Returns
    -------
    best_combination : dict
        dictionary of optimized transport values
    state : list
        list of temperature, pressure and mass_fractions corresponding to optimized transport values

    """

    species_names = mechanism.ctmech.species_names

    if not root:
        root = mechanism

    transport_variables = species_names + ['Pr']

    combination = {}

    var_error_square = []
    store_combination = []

    # Storing states
    store_temperatures = []
    store_pressures = []
    store_mass_fractions = []

    for case_for_data in cases_list:

        # Initializing the transport values with the burnt gases of first case
        if case_for_data.reactor_type != 'C1DP':
            temperatures = list(case_for_data.extract_quantity('T', mechanism=mechanism))
            temperature = max(temperatures)
            index_gas_state = temperatures.index(temperature)
            if type(index_gas_state) == list():
                index_gas_state = index_gas_state[0]
            indices_gas_state = [index_gas_state]
            list_of_states = ['burnt']
        else:
            indices_gas_state = [1, -1]
            list_of_states = ['fresh', 'burnt']

        for index_in_list, index_gas_state in enumerate(indices_gas_state):

            temperature = case_for_data.extract_quantity('T', mechanism=mechanism)[index_gas_state]
            store_temperatures.append(temperature)

            pressure = case_for_data.extract_quantity('P', mechanism=mechanism)[index_gas_state]
            store_pressures.append(pressure)

            mass_fractions = case_for_data.extract_quantity('Y', mechanism=mechanism)[index_gas_state]
            store_mass_fractions.append(mass_fractions)

            mechanism.ctmech.TPY = temperature, pressure, mass_fractions

            schmidt = tools.compute_Schmidt_numbers(mechanism.ctmech)
            prandtl = tools.compute_Prandtl_number(mechanism.ctmech)

            for index, var in enumerate(transport_variables):
                if var != 'Pr':
                    combination[var] = schmidt[index]
                else:
                    combination[var] = prandtl

            # Writing dummy mixture_database
            write_mixture_database(mechanism, transport_dict=combination)

            # Creating Mechanism object with AVBP transport
            mechanism_AVBP = mechanisms.Mechanism(mechanism.path, name=f'{mechanism.name}_AVBP_{case_for_data.myid}'
                                                     f'_{list_of_states[index_in_list]}', transport='AVBP')

            sum_square = 0
            for case_to_run in cases_list:
                # Loading reference case xml solution
                database.load_solutions_names(case_to_run, mechanism)
                file_to_restore = case_to_run.xml_solution

                case_to_run.run(mechanism_AVBP, restore=file_to_restore, overwrite=overwrite)

                error_list = error.case_error(case_to_run, root, mechanism=mechanism_AVBP)
                logger.info('')
                # Computing sum of the squared errors weighed by their tolerances
                weighed_error_list = [(error / tolerance) ** 2
                                      for error, tolerance in zip(error_list,
                                                                  list(case_to_run.error_dict.values()))]
                sum_square += sum(weighed_error_list)

            # Storing the errors
            var_error_square.append(sum_square)
            store_combination.append(combination.copy())

    index_best = var_error_square.index(min(var_error_square))

    # Setting the mechanism to the correct state
    state = [store_temperatures[index_best], store_pressures[index_best], store_mass_fractions[index_best]]
    mechanism.ctmech.TPY = store_temperatures[index_best], store_pressures[index_best], store_mass_fractions[index_best]

    best_combination = store_combination[index_best]

    return best_combination, state


def transport_error(cases_list, mechanism, best_combination, state, root=None,
                    output_dir="AVBP_input", save_init=True, plot=False):
    """Computes the error induced by AVBP transport

    Parameters
    ----------
    cases_list :
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    best_combination :
        dictionary of Schmidt and Prandtl values
    state :
        temperature, pressure and mass fractions corresponding to the best_combination set
    root :
        class :func:`~ARCANE.mechanisms.Mechanism` object to be compared to (Default value = None)
    output_dir :
        directory for storing the outputs (Default value = "AVBP_input")
    save_init :
        if True, saves initial composition for direct copy/paste in solutbound file (Default value = True)
    plot :
        if True, plots the error graph (Default value = False)


    Created: 18/12/15 [QC]

    """

    # Writing dummy mixture_database
    write_mixture_database(mechanism, transport_dict=best_combination, state=state)

    # Creating Mechanism object with AVBP transport
    mechanism_AVBP = mechanisms.Mechanism(mechanism.path, transport='AVBP')
    mechanism_AVBP.name += '_AVBP'

    if not root:
        root = mechanism

    for case in cases_list:
        # Loading reference case xml solution
        database.load_solutions_names(case, mechanism)
        file_to_restore = case.xml_solution

        case.run(mechanism_AVBP, restore=file_to_restore, overwrite=True)

        error.case_error(case, root, mechanism=mechanism_AVBP)
        logger.info('')

        if save_init:
            AVBP_init(case, mechanism_AVBP, sol_dir=output_dir + "/sol_can2av")
            write_inlet_for_solutbound(case, mechanism_AVBP, sol_dir=output_dir + "/solutbound_inlet")
            write_inlet_for_makesolution(case, mechanism_AVBP, sol_dir=output_dir + "/sol_makesolution")
            write_outlet_for_gas_out(case, mechanism_AVBP, sol_dir=output_dir + "/sol_gas_out")

    if plot:
        cases_list_per_fuel = {}

        mechanism.nickname = 'Mixture averaged transport'
        mechanism_AVBP.nickname = 'AVBP transport'

        for index, case in enumerate(cases_list):
            fuel = repr(case.fuel)
            if fuel not in cases_list_per_fuel:
                cases_list_per_fuel[fuel] = []

            if fuel in cases_list_per_fuel:
                cases_list_per_fuel[fuel].append(case)

        for index, fuel in enumerate(cases_list_per_fuel):
            ref_case = cases_list_per_fuel[fuel][0]
            graphs.parametric_plot(cases_list_per_fuel[fuel], [mechanism, mechanism_AVBP],
                                   'C1DP', 'Equivalence ratio', 'Sl',
                                   P=ref_case.pressure, T=ref_case.temperature, fuel=eval(fuel),
                                   hold=True, legend='Mixture averaged transport')


def write_mixture_database(mechanism, state=None, transport_dict=None, header='', mixture_name='', parsing=False,
                           fuel={}, oxidizer={}):
    """Writes the mixture_database.dat file for AVBP

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    state :
        list of form [temperature, pressure, mass_fractions] to define the gas state (Default value = None)
    transport_dict :
        user-defined directory with species Schmidt numbers and Prandtl number (Default value = None)
    header :
        mixture_database.dat header (Default value = '')
    mixture_name :
        name of the mixture (Default value = '')
    parsing :
        if True, adds keywords for AVBP parsing (Default value = False)
    fuel :
        fuel composition (dictionary if several molecules) (Default value = {})
    oxidizer :
        oxydizer composition (dictionary if several molecules) (Default value = {})


    Created: 18/12/15 [QC]

    """

    ctmech = mechanism.ctmech
    if state:
        ctmech.TPY = state[0], state[1], state[2]

    # Initialization of the Transport Class
    ctmech.transport_model = 'Mix'
    if not mixture_name:
        mixture_name = ctmech.name

    # Reference values
    ref_temperature = format_to_AVBP(ctmech.T)
    ref_viscosity = format_to_AVBP(ctmech.viscosity)

    if transport_dict:
        Sch = []
        for spec in ctmech.species_names:
            Sch.append(transport_dict[spec])
        Pr = transport_dict['Pr']
    else:
        # Computing Schmidt and Prandtl
        Sch = tools.compute_Schmidt_numbers(ctmech)
        Pr = tools.compute_Prandtl_number(ctmech)

    # Converting lists and values to correct strings
    Sch_string = ''
    for index in range(ctmech.n_species):
        Sch_string += ' ' + format_to_AVBP(Sch[index])

    Pr_string = format_to_AVBP(Pr)

    species_string = ''
    for spec in ctmech.species_names:
        spec_string = ' ' + spec
        species_string += spec_string

    viscosity_exponent = format_to_AVBP(compute_power_law_exponent(mechanism, 'Mix'))

    if parsing:
        parsing_start = "$MIXTURE\n"
        parsing_end = "$end_MIXTURE"
    else:
        parsing_start = ""
        parsing_end = ""

    if type(fuel) is str:
        fuel = {fuel: 1}

    if type(oxidizer) is str:
        oxidizer = {oxidizer: 1}

    fuel_string = ''
    for spec in ctmech.species_names:
        if spec in fuel:
            fuel_string += ' ' + format_to_AVBP(fuel[spec])
        else:
            fuel_string += ' ' + format_to_AVBP(0.)

    oxi_string = ''
    for spec in ctmech.species_names:
        if spec in oxidizer:
            oxi_string += ' ' + format_to_AVBP(oxidizer[spec])
        else:
            oxi_string += ' ' + format_to_AVBP(0.)

        template = f"""\
{parsing_start}
!--------------------------------------------------------------------------------------------------
! {header}
!--------------------------------------------------------------------------------------------------
mixture_name = {mixture_name}
transport = computed
combustion_chemistry = analytic
species_name = {species_string}
species_Schmidt_number =  {Sch_string}
viscosity_law = power
mu_ref =   {ref_viscosity}
T_ref =   {ref_temperature}
viscosity_law_coeff =  {viscosity_exponent}
prandtl_number =   {Pr_string}
fuel_tank_composition = {fuel_string}
oxydizer_tank_composition = {oxi_string}

{parsing_end}
"""

    f = open('mixture_database.dat', 'w+')
    f.write(template)
    f.close()

    logger.info('mixture_database.dat written successfully')


def write_species_database(mechanism, fuel=[], liquid=[], liquid_properties={}, max_temperature=5000):
    """Writes the species_database.dat file for AVBP

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    fuel :
        list of fuel species (Default value = [])
    liquid :
        list of liquid species (Default value = [])
    liquid_properties :
        dictionary with liquid properties in the form
        {'species': {'liquid_density': 0,
        'boiling_temperature': 0
        'critical_temperature': 0
        'critical_pressure': 0
        'critical_volume': 0
        'viscosity_params': [0, 0, 0, 0]}} (Default value = {})
    max_temperature :
        maximum temperature of the table (Default value = 5000)


    Created: 18/12/15 [QC]
    Modified: 2O/03/24 [JW]

    """

    # Reference enthalpy of atoms at 0K
    atomic_h = {'C': -1054, 'H': -4234, 'O': -4340, 'N': -4335}

    ctmech = mechanism.ctmech
    species = ctmech.species()
    elements = ctmech.element_names
    mol_weights = ctmech.molecular_weights
    formation_enthalpies = []
    for spec in species:
        spec_enthalpy = spec.thermo.h(1) / 1000
        # Correcting the enthalpy as Cantera is set for 298.15 K
        for element in elements:
            if element in atomic_h:
                spec_enthalpy -= ctmech.n_atoms(spec.name, element) * atomic_h[element]
            else:
                if ctmech.n_atoms(spec.name, element) > 0:
                    logger.error(
                        f'The enthalpy of formation conversion from 298.15 K to 0K could not be done correctly.\n'
                        f'Atomic value of {element} for species {spec} is not known.'
                        f'Modifications in the code are required (AVBP.write_species_database function).')

        if abs(spec_enthalpy) < 1e3:
            spec_enthalpy = 0

        formation_enthalpies.append(spec_enthalpy)

    ctmech.TP = 300, 1e5
    ctmech.transport_model = 'Mix'
    species_viscosities = ctmech.species_viscosities

    f = open('species_database.dat', 'w+')

    f.write('$SPECIES ! Gaseous, solid and liquid species')

    template = Template("""\

!---------------------------
!
!    $name
!
!---------------------------
species_name = $name
species_is_fuel = $isfuel
species_type = $state
species_C_atoms = $carbon_number
species_H_atoms = $hydrogen_number
species_O_atoms = $oxygen_number
species_N_atoms = $nitrogen_number
species_molecular_weight =   $molecular_weight
species_formation_enthalpy =   $formation_enthalpy
species_ref_molecular_viscosity =   $ref_mol_viscosity
species_ref_temperature =   $ref_temperature
species_power_law_exponent =   $power_law_exponent
species_enthalpy_table =  $enthalpy_table
species_entropy_table =   $entropy_table
""")

    template_tpf = Template("""\
tpf_species_density = $tpf_density
tpf_species_surface_tension = $tpf_surface_tension
tpf_species_Tcrit = $tpf_critical_temperature
tpf_species_Vcrit = $tpf_critical_volume
tpf_species_Hlat = $tpf_latent_heat
tpf_species_Psat = $tpf_saturation_pressure
tpf_species_viscosity = $tpf_viscosity
""")

    for index, spec in enumerate(species):

        name = spec.name
        if name in fuel:
            isfuel = 'yes'
        else:
            isfuel = 'no'
        if name in liquid:
            state = 'gaseous_and_liquid'
        else:
            state = 'gaseous'

        if 'C' in ctmech.element_names:
            carbon_number = ctmech.n_atoms(spec.name, 'C')
        else:
            carbon_number = 0

        if 'H' in ctmech.element_names:
            hydrogen_number = ctmech.n_atoms(spec.name, 'H')
        else:
            hydrogen_number = 0

        if 'O' in ctmech.element_names:
            oxygen_number = ctmech.n_atoms(spec.name, 'O')
        else:
            oxygen_number = 0

        if 'N' in ctmech.element_names:
            nitrogen_number = ctmech.n_atoms(spec.name, 'N')
        else:
            nitrogen_number = 0

        molecular_weight = mol_weights[index] / 1000

        formation_enthalpy = formation_enthalpies[index]

        ref_mol_viscosity = species_viscosities[index]

        ref_temperature = 300.0

        power_law_exponent = compute_power_law_exponent(mechanism, spec.name)

        enthalpy_table = generate_thermo_table(mechanism, 'enthalpy', spec.name, max_temperature=max_temperature)

        entropy_table = generate_thermo_table(mechanism, 'entropy', spec.name, max_temperature=max_temperature)

        species_string = template.substitute(name=name,
                                             isfuel=isfuel,
                                             state=state,
                                             carbon_number=format_to_AVBP(carbon_number),
                                             hydrogen_number=format_to_AVBP(hydrogen_number),
                                             oxygen_number=format_to_AVBP(oxygen_number),
                                             nitrogen_number=format_to_AVBP(nitrogen_number),
                                             molecular_weight=format_to_AVBP(molecular_weight),
                                             formation_enthalpy=format_to_AVBP(formation_enthalpy),
                                             ref_mol_viscosity=format_to_AVBP(ref_mol_viscosity),
                                             ref_temperature=format_to_AVBP(ref_temperature),
                                             power_law_exponent=format_to_AVBP(power_law_exponent),
                                             enthalpy_table=format_to_AVBP_list(enthalpy_table),
                                             entropy_table=format_to_AVBP_list(entropy_table))

        tpf_density = 0
        tpf_surface_tension = 0
        tpf_critical_temperature = 0
        tpf_critical_volume = 0
        tpf_latent_heat = 0
        tpf_saturation_pressure = 0
        tpf_viscosity = 0

        # If they exist, retrieve the liquid properties
        if name in liquid_properties:
            # From input
            species_tpf_dict = liquid_properties[name]
        elif name in liquid:
            # From stored database
            species_tpf_dict = chem_data.get_tpf_data(name)
        else:
            species_tpf_dict = {}

        if species_tpf_dict:
            species_tpf_dict['molecular_weight'] = molecular_weight

            tpf_density = species_tpf_dict['liquid_density']
            tpf_critical_temperature = species_tpf_dict['critical_temperature']
            tpf_critical_volume = species_tpf_dict['critical_volume']

            # Computing interesting parameters
            # Every expression is from
            # Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).
            # The properties of gases and liquids. New York: McGraw-Hill.
            tpf_saturation_pressure = species_sat_pressure(species_tpf_dict)
            tpf_latent_heat = species_latheat_vapor(species_tpf_dict)
            tpf_viscosity = species_liquid_viscosity(species_tpf_dict)
            tpf_surface_tension = species_surface_tension(species_tpf_dict)

        species_string_tpf = template_tpf.substitute(tpf_density=format_to_AVBP(tpf_density),
                                                     tpf_critical_temperature=format_to_AVBP(tpf_critical_temperature),
                                                     tpf_critical_volume=format_to_AVBP(tpf_critical_volume),
                                                     tpf_saturation_pressure=format_to_AVBP(tpf_saturation_pressure),
                                                     tpf_latent_heat=format_to_AVBP(tpf_latent_heat),
                                                     tpf_viscosity=format_to_AVBP(tpf_viscosity),
                                                     tpf_surface_tension=format_to_AVBP(tpf_surface_tension))

        f.write(species_string)
        f.write(species_string_tpf)
        f.write('!---------------------------\n')

    f.write('$end_SPECIES')
    f.close()

    logger.info('species_database.dat written successfully')


def format_to_AVBP(value, decimals=8):
    """Takes a float and convert it with the right AVBP formatting

    Parameters
    ----------
    value :
        float to convert
    decimals :
        number of significant numbers (Default value = 8)

    Returns
    -------
    string : str
        formatted string


    Created: 18/12/15 [QC]

    """

    if isinstance(value, list):
        string = format_to_AVBP_list(value, decimals=decimals)
    else:
        string = format(value, '.' + str(decimals) + 'e')
        string = string.replace('e+', 'd')
        string = string.replace('e-', 'd-')

    return string


def format_to_AVBP_list(values_list, decimals=8):
    """Takes a list of float and convert it with the right AVBP formatting

    Parameters
    ----------
    values_list :
        list of floats to convert
    decimals :
        number of significant numbers (Default value = 8)

    Returns
    -------
    string_list : str
        formatted string


    Created: 18/12/15 [QC]

    """

    string_list = ''
    for value in values_list:
        string = format_to_AVBP(value, decimals=decimals)

        string_list += '  ' + string

    return string_list


def generate_thermo_table(mechanism, quantity, species, max_temperature=5000):
    """Generates thermodynamic tables

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    quantity :
        entropy or enthalpy
    species :
        name of concerned species
    max_temperature :
        maximum temperature of the table (Default value = 5000)

    Returns
    -------
    table : list
        thermodynamic table


    Created: 18/12/15 [QC]

    """

    ctmech = mechanism.ctmech
    table = []
    temperatures = list(range(0, max_temperature + 100, 100))
    temperatures[0] = 1

    # Reference temperature to 1 K as 0 K is not possible
    ctmech.TP = 1, 1e5
    correction_enthalpy = ctmech.partial_molar_enthalpies
    correction_entropy = ctmech.partial_molar_entropies

    spec_index = ctmech.species_index(species)

    sliding_array = [0, -1, -0.2]
    linear = False
    constant_increment = 0
    for index, temperature in enumerate(temperatures):

        ctmech.TP = temperature, 1e5
        if not linear:
            if quantity == 'enthalpy':
                new_enthalpy = ctmech.partial_molar_enthalpies[spec_index] - correction_enthalpy[spec_index]
                table.append(new_enthalpy / 1000)
            elif quantity == 'entropy':
                new_entropy = ctmech.partial_molar_entropies[spec_index] - correction_entropy[spec_index]
                table.append(new_entropy / 1000)
            else:
                logger.warning('Thermodynamic quantity not valid (enthalpy or entropy)')
                break

            # Correction of decreasing enthalpy by linear prolongation
            sliding_array[0] = sliding_array[1]
            sliding_array[1] = sliding_array[2]
            sliding_array[2] = table[-1]

            slope = np.gradient(sliding_array, [0, 100, 200])
            curve = np.gradient(slope, [0, 100, 200])

            # Sets a linear continuity is curve is concave
            # And temperature higher than lower threshold of NASA polynomial
            if curve[1] <= 0 and temperature - 200 > 200:
                linear = True
                constant_increment = slope[0] * 100
        else:
            previous_value = table[-1]
            table.append(previous_value + constant_increment)

    return table


def compute_power_law_exponent(mechanism, mixture):
    """Computes the exponent of the viscosity power law

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    mixture :
        mix for the whole mixture or species name

    Returns
    -------
    exponent : float
        power law exponent


    Created: 18/12/15 [QC]

    """

    ctmech = mechanism.ctmech

    T_ref = ctmech.T
    P = ctmech.P
    Y = ctmech.Y

    ctmech.transport_model = 'Mix'
    species_names = ctmech.species_names

    if mixture in species_names:
        mu_ref = ctmech.species_viscosities[ctmech.species_index(mixture)]
    elif mixture in ['Mix', 'mix']:
        mu_ref = ctmech.viscosity

    # Set constant parameters and shortcuts for the computation of viscosity
    Tstart = 200
    Tend = 10000
    n_temperatures = 100

    # Temperature array, equidistant in T
    temperatures = []
    mu = []
    mu_norm = []
    temperatures_norm = []

    dx = (Tend - Tstart) / n_temperatures
    for i in range(n_temperatures):
        # Temperatures list
        new_T = Tstart + dx * i
        temperatures.append(new_T)

        if mixture in species_names:
            mu_new = ctmech.species_viscosities[ctmech.species_index(mixture)]
        elif mixture in ['Mix', 'mix']:
            mu_new = ctmech.viscosity

        # Exact viscosity list
        ctmech.TPY = new_T, P, Y
        mu.append(mu_new)

        # mu_norm = mu / mu_ref
        mu_norm.append(mu_new / mu_ref)
        # T_norm = T / T_ref
        temperatures_norm.append(ctmech.T / T_ref)

    # Fitting
    exponent, pcov = opt.curve_fit(power, temperatures_norm, mu_norm)
    exponent = float(exponent)

    return exponent


def power(x, n):
    """Computes x^n

    Parameters
    ----------
    x :
        variable
    n :
        exponent

    Returns
    -------
    x ** n : float
        x^n

    """
    return x ** n


def AVBP_init(case, mechanism, sol_dir="sol_can2av"):
    """Writes the csv for 1D flame initialization in AVBP

    Parameters
    ----------
    case :
        class :func:`~ARCANE.cases.Case` object
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    sol_dir :
         directory where to store solutions (Default value = "sol_can2av")


    Created: 19/01/25 [QC]

    """

    # Create it if no already there (do not overwrite)
    database.create_dir(sol_dir, overwrite=False)

    # Compiling mechanism (does nothing if mechanism is not ARC)
    mechanism.compile()

    f = ct.FreeFlame(mechanism.ctmech)

    # Loading reference case xml solution
    database.load_solutions_names(case, mechanism)
    file_to_restore = case.xml_solution

    solution_path, solution_name = database.xml_solution_file_and_name(file_to_restore)
    f.restore(solution_path, name=solution_name, loglevel=0)
    file_name = sol_dir + "/" + case.myid + ".csv"

    # Impose sum of mass fractions equal to 1
    logger.info('Normalization to impose sum of mass fractions equal to 1')
    sumYk = np.sum(f.Y, 0)
    for k, species in enumerate(f.gas.species_names):
        f.set_profile(species, f.grid / f.grid[-1], f.Y[k] / sumYk)

    try:
        f.write_AVBP(file_name)
    except BrokenPipeError:
        logger.warning('WARNING !!!! Failed to write csv file as ' + file_name)
    else:
        logger.debug('Cantera initial solution for AVBP written as ' + file_name)

    # Freeing dynamic library
    mechanism.reset()

    # write can2av files
    cantera_sol = './' + case.myid + ".csv"
    AVBP_sol = './' + case.myid + '_AVBP.h5'
    name_can2av = sol_dir + "/" + "can2av.choices_" + case.myid

    file = open(name_can2av, "w")
    file.write(f"'{cantera_sol}'            ! CANTERA last or restart file\n")
    file.write(f"'{AVBP_sol}'          ! The AVBP solution\n")
    file.write("'../MESH/mesh.mesh.h5'                ! The AVBP mesh file\n")
    file.write(f"{f.flame.P}                        ! Pressure of CANTERA solution\n")
    file.write(f"{mechanism.ctmech.n_species}                     ! Number of species in CANTERA solution\n")
    file.write(f"{len(f.flame.grid)}                              ! Number of grid points in CANTERA solution\n")
    file.write(f"{0}                                ! Number of reactions in CANTERA solution")
    file.close()


def write_inlet_for_solutbound(case, mechanism, sol_dir="solutbound_inlet"):
    """Writes initial composition for direct copy/paste in solutbound file

    Parameters
    ----------
    case :
        class :func:`~ARCANE.cases.Case` object
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    sol_dir :
        directory in which files will be stored (Default value = "solutbound_inlet")


    Created: 19/01/25 [QC]

    """

    file_name = case.myid + "_inlet.txt"

    # Create it if no already there (do not overwrite)
    database.create_dir(sol_dir, overwrite=False)

    data_dict = case.data_dict(mechanism)

    sum = 0
    mass_fractions_list = []

    f = open(sol_dir + "/" + file_name, 'w+')

    for spec in mechanism.ctmech.species_names:

        mass_fraction = round(float(data_dict[spec][0]), 5)
        if mass_fraction > 1e-15:
            mass_fractions_list.append(mass_fraction)
        else:
            mass_fractions_list.append(0)

        text = spec + "\n" + format_to_AVBP(mass_fraction, decimals=4) + "\n\n"

        sum += mass_fraction

        f.write(text)

    f.close()

    if sum != 1:

        f = open(sol_dir + "/" + file_name, 'w+')

        for index, spec in enumerate(mechanism.ctmech.species_names):
            mass_fractions_list[index] /= sum

            text = spec + "\n" + format_to_AVBP(mass_fractions_list[index], decimals=4) + "\n\n"

            f.write(text)

        f.close()


def write_outlet_for_gas_out(case, mechanism, sol_dir="sol_gas_out"):
    """Writes burnt gases composition for copy/past in gas_out.choices

    Parameters
    ----------
    case :
        class :func:`~ARCANE.cases.Case` object
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    sol_dir :
        directory in which files will be stored (Default value = "sol_gas_out")


    Created: 19/01/25 [QC]

    """

    file_name = case.myid + "_gas_out.txt"

    # Create it if no already there (do not overwrite)
    database.create_dir(sol_dir, overwrite=False)

    data_dict = case.data_dict(mechanism)

    adiabatic_temperature = data_dict['T'][-1]

    sum = 0
    mass_fractions_list = []

    f = open(sol_dir + "/" + file_name, 'w+')

    f.write(format_to_AVBP(adiabatic_temperature, decimals=4) + '                          ! T_burnt\n')

    for spec in mechanism.ctmech.species_names:

        mass_fraction = round(float(data_dict[spec][-1]), 5)
        if mass_fraction > 1e-15:
            mass_fractions_list.append(mass_fraction)
        else:
            mass_fractions_list.append(0)

        text = "'" + spec + "' " + format_to_AVBP(mass_fraction, decimals=4) + " 0 0"

        if spec == mechanism.ctmech.species_names[0]:
            text += '                ! Species, mass fractions in burnt gases, fuel tank and oxy. tank\n'
        else:
            text += '\n'

        sum += mass_fraction

        f.write(text)

    f.close()

    if sum != 1:

        f = open(sol_dir + "/" + file_name, 'w+')

        f.write(format_to_AVBP(adiabatic_temperature, decimals=4) + '                           ! T_burnt\n')

        for index, spec in enumerate(mechanism.ctmech.species_names):

            mass_fractions_list[index] /= sum

            text = "'" + spec + "' " + format_to_AVBP(mass_fractions_list[index], decimals=4) + " 0 0"

            if index == 0:
                text += '                ! Species, mass fractions in burnt gases, fuel tank and oxy. tank\n'
            else:
                text += '\n'

            f.write(text)

        f.close()


def write_inlet_for_makesolution(case, mechanism, sol_dir="sol_makesolution"):
    """Writes fresh gases composition for copy/past in makesolution.choices

    Parameters
    ----------
    case :
        class :func:`~ARCANE.cases.Case` object
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    sol_dir :
        directory in which files will be stored (Default value = "sol_makesolution")


    Created: 19/01/25 [QC]

    """

    file_name = case.myid + "_makesolution.txt"

    # Create it if no already there (do not overwrite)
    database.create_dir(sol_dir, overwrite=False)

    data_dict = case.data_dict(mechanism)

    pressure = data_dict['P'][0]
    temperature = data_dict['T'][0]

    sum = 0

    f = open(sol_dir + "/" + file_name, 'w+')

    f.write(format_to_AVBP(pressure, decimals=4) + '          ! P\n')
    f.write(format_to_AVBP(temperature, decimals=4) + '          ! T\n')
    f.write('0.0d0              ! U\n')
    f.write('0.0d0              ! V\n')
    f.write('0.0d0              ! W\n')

    mass_fractions_list = []
    for spec in mechanism.ctmech.species_names:

        mass_fraction = round(float(data_dict[spec][0]), 5)
        if mass_fraction > 1e-15:
            mass_fractions_list.append(mass_fraction)
        else:
            mass_fractions_list.append(0)

        text = format_to_AVBP(mass_fraction, decimals=4) + "          ! " + spec + '\n'

        sum += mass_fraction

        f.write(text)

    f.close()

    if sum != 1:

        f = open(sol_dir + "/" + file_name, 'w+')

        f.write(format_to_AVBP(pressure, decimals=4) + '          ! P\n')
        f.write(format_to_AVBP(temperature, decimals=4) + '          ! T\n')
        f.write('0.0d0            ! U\n')
        f.write('0.0d0            ! V\n')
        f.write('0.0d0            ! W\n')

        for index, spec in enumerate(mechanism.ctmech.species_names):
            mass_fractions_list[index] /= sum

            text = format_to_AVBP(mass_fractions_list[index], decimals=4) + "          ! " + spec + '\n'

            f.write(text)

        f.close()


def temperature_array_parameters(species_tpf_dict):
    """Constructs a dictionary with the correct temperature array to fit data on
    The array corresponds to one value each 10K from 0K to the critical temperature

    Parameters
    ----------
    species_tpf_dict :
        dictionary with liquid properties in the form
        {'liquid_density': 0,
        'boiling_temperature': 0
        'critical_temperature': 0
        'critical_pressure': 0
        'critical_volume': 0
        'viscosity_params': [0, 0, 0, 0]}

    Returns
    -------
    parameters_dict : dict
        dictionary containing the parameters (temperature array, temperature ratio, size of the array

    """

    critical_temperature = species_tpf_dict['critical_temperature']
    boiling_temperature = species_tpf_dict['boiling_temperature']

    # Declaring size of temperature fit (from 0K to critical temperature)
    array_size = int(critical_temperature // 10) + int(2)
    temperature_array = np.linspace(0, 10 * (array_size - 1), num=array_size)
    temperature_ratio_array = temperature_array / critical_temperature
    boiling_critical_ratio = boiling_temperature / critical_temperature

    parameters_dict = {'temperature_array': temperature_array,
                       'temperature_ratio_array': temperature_ratio_array,
                       'boiling_critical_ratio': boiling_critical_ratio,
                       'array_size': array_size}

    return parameters_dict


def species_sat_pressure(species_tpf_dict):
    """Calculates the saturation pressure of a species

    Parameters
    ----------
    species_tpf_dict :
        dictionary with liquid properties in the form
        {'liquid_density': 0,
        'boiling_temperature': 0
        'critical_temperature': 0
        'critical_pressure': 0
        'critical_volume': 0
        'viscosity_params': [0, 0, 0, 0]}

    Returns
    -------
    psat_out : list
        saturation pressure array

    """

    temperature_parameters = temperature_array_parameters(species_tpf_dict)

    array_size = temperature_parameters['array_size']
    temperature_ratio_array = temperature_parameters['temperature_ratio_array']
    boiling_critical_ratio = temperature_parameters['boiling_critical_ratio']

    critical_pressure = species_tpf_dict['critical_pressure']

    psat_calc = np.zeros(array_size)

    # Interesting variables for the computing
    tau_array = 1 - temperature_ratio_array
    taub = 1 - boiling_critical_ratio

    tau_array = tau_array[1:array_size - 1]
    temperature_ratio_array = temperature_ratio_array[1:array_size - 1]

    # Calculation of the eccentric factor
    f0b = (-5.97616 * taub + 1.29874 * (taub ** 1.5)
           - 0.60394 * (taub ** 2.5) - 1.06841 * (taub ** 5)) / boiling_critical_ratio
    f1b = (-5.03365 * taub + 1.11505 * (taub ** 1.5)
           - 5.41217 * (taub ** 2.5) - 7.46628 * (taub ** 5)) / boiling_critical_ratio
    omega = - (np.log(critical_pressure / 1.01325) + f0b) / f1b

    # Calculation of the saturation pressure
    f0 = (-5.97616 * tau_array + 1.29874 * (tau_array ** 1.5)
          - 0.60394 * (tau_array ** 2.5) - 1.06841 * (tau_array ** 5)) / temperature_ratio_array
    f1 = (-5.03365 * tau_array + 1.11505 * tau_array ** 1.5
          - 5.41217 * tau_array ** 2.5 - 7.46628 * tau_array ** 5) / temperature_ratio_array
    f2 = (-0.64771 * tau_array + 2.41539 * tau_array ** 1.5
          - 4.26979 * tau_array ** 2.5 + 3.25259 * tau_array ** 5) / temperature_ratio_array

    psat_calc[1:array_size - 1] = np.exp(f0 + f1 * omega + f2 * (omega ** 2)) * critical_pressure * 1e5

    # Output formatting
    psat_out = psat_calc.tolist()

    return psat_out


def species_latheat_vapor(species_tpf_dict):
    """Calculates the latent heat of vaporization of a species

    Parameters
    ----------
    species_tpf_dict :
        dictionary with liquid properties in the form
        {'liquid_density': 0,
        'boiling_temperature': 0
        'critical_temperature': 0
        'critical_pressure': 0
        'critical_volume': 0
        'viscosity_params': [0, 0, 0, 0]}

    Returns
    -------
    hlat_out : list
        latent heat of vaporization array

    """

    temperature_parameters = temperature_array_parameters(species_tpf_dict)

    array_size = temperature_parameters['array_size']
    temperature_ratio_array = temperature_parameters['temperature_ratio_array']
    boiling_critical_ratio = temperature_parameters['boiling_critical_ratio']

    critical_pressure = species_tpf_dict['critical_pressure']
    boiling_temperature = species_tpf_dict['boiling_temperature']
    molecular_weight = species_tpf_dict['molecular_weight']

    hlat_calc = np.zeros(array_size)

    # Calculation of the latent heat of vaporization
    hlat_tb = 8.314 * boiling_temperature \
              * ((1 - boiling_critical_ratio) ** 0.38) * (np.log(critical_pressure) - 0.513 + 0.5066
                                                          / (
                                                                  critical_pressure * boiling_critical_ratio ** 2)) / molecular_weight

    hlat_tb = hlat_tb / (1 - boiling_critical_ratio
                         + (1 - (1 - boiling_critical_ratio) ** 0.38) * np.log(boiling_critical_ratio))

    scale_const = hlat_tb / (1 - boiling_critical_ratio) ** 0.38

    hlat_calc[0:array_size - 1] = scale_const * (1 - temperature_ratio_array[0:array_size - 1]) ** 0.38

    # Output formatting
    hlat_out = hlat_calc.tolist()

    return hlat_out


def species_liquid_viscosity(species_tpf_dict):
    """Calculates the liquid viscosity of a species

    Parameters
    ----------
    species_tpf_dict :
        dictionary with liquid properties in the form
        {'liquid_density': 0,
        'boiling_temperature': 0
        'critical_temperature': 0
        'critical_pressure': 0
        'critical_volume': 0
        'viscosity_params': [0, 0, 0, 0]}

    Returns
    -------
    eta_out : list
        liquid viscosity array

    """

    temperature_parameters = temperature_array_parameters(species_tpf_dict)

    array_size = temperature_parameters['array_size']
    temperature_array = temperature_parameters['temperature_array']

    viscosity_params = species_tpf_dict['viscosity_params']

    eta_calc = np.zeros(array_size)

    # Viscosity calculation
    eta_calc[20:array_size - 1] = 0.001 * 10 ** (
            viscosity_params[0] + viscosity_params[1] / temperature_array[20:array_size - 1] +
            viscosity_params[2] * temperature_array[20:array_size - 1] +
            viscosity_params[3] * temperature_array[20:array_size - 1] ** 2.0)

    # Output formatting
    eta_out = eta_calc.tolist()

    return eta_out


def species_surface_tension(species_tpf_dict, calculation_temperature=350):
    """
    Calculates the surface tension of a species

    Parameters
    ----------
    species_tpf_dict :
        dictionary with liquid properties in the form
        {'liquid_density': 0,
        'boiling_temperature': 0
        'critical_temperature': 0
        'critical_pressure': 0
        'critical_volume': 0
        'viscosity_params': [0, 0, 0, 0]}
    calculation_temperature: float
        temperature for the calculation

    Returns
    -------
    sigma_out : list
        surface tension array

    """
    calculation_temperature = int(calculation_temperature / 10)
    temperature_parameters = temperature_array_parameters(species_tpf_dict)

    temperature_ratio_array = temperature_parameters['temperature_ratio_array']
    boiling_critical_ratio = temperature_parameters['boiling_critical_ratio']

    critical_pressure = species_tpf_dict['critical_pressure']
    critical_temperature = species_tpf_dict['critical_temperature']

    # Calculation of the surface tension at 350K (35th value of the liquid temperature vector)
    q_visc = 0.1196 * (1 + (boiling_critical_ratio * np.log(critical_pressure / 1.01325)) / (
            1 - boiling_critical_ratio)) - 0.279

    sigma_calc = (critical_pressure ** (2 / 3)) * (critical_temperature ** (1 / 3)) * q_visc * (
            (1 - temperature_ratio_array[calculation_temperature]) ** (11 / 9)) * 0.001

    # Output formatting
    sigma_out = sigma_calc

    return sigma_out
