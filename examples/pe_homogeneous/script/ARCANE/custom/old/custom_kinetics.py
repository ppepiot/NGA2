"""
Utilities to write fortran source codes for a variety of outlets: Cantera, AVBP, YALES2 etc.
Includes QSS handling
"""

import filecmp
import os.path
import shutil
import subprocess
import sys
import time

import ARCANE.ct2cti as ct2cti
import ARCANE.custom.custom_printing as pr
import ARCANE.custom.custom_qss as qss
import ARCANE.display as display
import ARCANE.database as database
import ARCANE.mechanisms as mechanisms
import ARCANE.networks as networks
import ARCANE.tools as tools
import cantera as ct
import numpy as np
import scipy.optimize as opt

logger = display.Logger()
logger.set_log('logCompute')


def print_fortran(mechanism, f90_filename='custom_kinetics.f90', species_qss_names=[], use='Cantera',
                  routine_name='customkinetics',semi_implicit=False, rrate=True, force=True,
                  exponential=False,implicit_exponential=False,classical=True,dt_exp_user=1e-10,
                  supplement_dict=None, maintain_given_order=True, fit_order=0,
                  start_temperature=300, end_temperature=2500, number_of_temperatures=6000, pressure=101325,
                  verbose=False):

    """Comprehensive function to output fortran source code, called from outside.

    Parameters
    ----------
    mechanism :
        `ARCANE.mechanisms.Mechanism` object that contains all mechanism-specific info
    f90_filename : str
        name of the output .f90 file (Default value = 'custom_kinetics.f90')
    species_qss_names : list
        list of the QSS species. (Default value = [])
    use : str
        code that will be using the fortran routine: 'Cantera' (Default), 'AVBP', 'YALES2', 'PASR'
    routine_name : str
        name of the fortran routine inside fortran file (Default value = 'customkinetics')
    semi_implicit : list
        list of species for writing the semi-implicitation routine for AVBP (Default value = False)
    rrate : bool
        if True, writes the rrate routine for AVBP (Default value = True)
    force : bool
        if True, species and reactions will be removed from the mechanism if they induce QSS error (Default value = True)
    exponential : bool
        if True uses the exponential formulation for production rates (only for AVBP) (Default value = False)
    implicit_exponential :
        if True uses the implicit exponential formulation for production rates (only for `NTMIX` now) (Default value = False)
    classical : bool
        if True add routine for classical formulation of production rates (only for `NTMIX` now) (Default value = True)
    dt_exp_user : flat
        only used if exponential=True, and for Cantera. The value is the fixed dt used for exponential integration. (Default value = 1e-10)
    supplement_dict : dict
        dictionary containing the infos about the mechanism that will be written in the header (Default value = None):
        - `'details'`: information about applications (fuel/oxidizer)
        - `'authors'`: name of the persons generating the file
        - `'since'`: version of the CFD code used
        - `'note'`: everything else (range of validity, specificity, etc.)
    maintain_given_order : bool
        if True, the species inside the groups will be sorted according to the global list (Default value = True)
    fit_order : int
        order of the polynomial used for fitting the reverse arrhenius functions (Default value = 0)
    start_temperature : float
        start temperature for the fit (Default value = 300)
    end_temperature : float
        end temperature for the fit (Default value = 2500)
    number_of_temperatures : int
        number of temperatures for the fit (Default value = 6000)
    pressure : float
        pressure used in the equilibrium constant computation (Default value = 101325)
    verbose : bool
        if True, writes a message when file is written (Default value = False)

    """
    # Checks options compatibility
    if use == "NTMIX":
        if exponential and implicit_exponential:
            logger.error('exponential and implicit_exponential options cannot be used at the same time.')
            logger.error('Please deactivate one of the options. ')
            sys.exit()
        if not exponential and not classical:
            classical = True
            logger.error('Classical chemistry used by default')
        if implicit_exponential:
            classical = True
            logger.error('Classical chemistry used by default with implicit exponential')

    # Ideal gas constant, converted to SI from Cantera's J/kmol/K
    Rcst = ct.gas_constant / 1000.0

    # Checks if the mechanism is an arc or a skeletal
    # And loads the Cantera Solution
    if mechanism.nqss > 0:
        species_qss_names = mechanism.species_qss_names
        mechanism = mechanism.skeletal

    ctmech = mechanism.ctmech

    # Gets lists of species and reactions
    elements = ctmech.element_names

    species = ctmech.species
    species_names = ctmech.species_names
    reactions = ctmech.reactions()

    n_atoms_in_spec = ctmech.n_atoms

    ne = len(elements)

    # Checking for inert species
    nu = abs(mechanism.network.nup - mechanism.network.nur)
    inert = [spec for index, spec in enumerate(species_names) if sum(nu[index, :]) == 0]

    mono_element = []
    mono_element_species = []

    for spec in inert:
        spec_obj = species(spec)
        if len(spec_obj.composition) == 1:
            mono_element.append(list(spec_obj.composition.keys())[0])
            mono_element_species.append(spec)

    ne -= len(mono_element)
    [elements.remove(ele) for ele in mono_element if ele in elements]

    ns = len(species_names)
    nr = len(reactions)
    timescales = mechanism.timescales

    # =================================================================
    # Initialization - QSS species treatment
    # =================================================================

    if force:
        remove_from = 'mechanism'
    else:
        remove_from = 'list'

    # index of last species (nqss-1)
    nqss = len(species_qss_names)

    # Skip if no qss
    if nqss > 0:

        # Check that all species in QSS list are indeed valid species in mechanism
        qss_not_found = [spec for spec in species_qss_names if spec not in species_names]
        if qss_not_found:
            for qss_spec in qss_not_found:
                logger.error('Species ' + qss_spec + ' not in mechanism')
            quit()

        new_ctmech, new_species_qss_names = problematic_qss(mechanism, species_qss_names,
                                                            remove_from=remove_from)

        if new_species_qss_names != species_qss_names:
            species_qss_names = new_species_qss_names
            mechanism = mechanisms.Mechanism(new_ctmech, parent=mechanism.parent, how=mechanism.how)
            mechanism.timescales = timescales

        # Reload mechanism in Cantera, and get lists of species and reactions
        ctmech = mechanism.ctmech
        species_names = ctmech.species_names
        reactions = ctmech.reactions()
        ns = len(species_names)
        nr = len(reactions)
        nqss = len(species_qss_names)

    # Correct naming for fortran while saving a copy congruent with Cantera
    species_names_original = species_names.copy()
    species_names = tools.convert_to_valid_fortran_var(species_names)
    species_qss_names_original = species_qss_names.copy()
    species_qss_names = tools.convert_to_valid_fortran_var(species_qss_names)

    # Establish handle to fortran source file
    if use == 'NTMIX':
        if not f90_filename.endswith('.F90'):
            if f90_filename.endswith('.f90'):
                f90_filename = f90_filename[:-4]

            f90_filename += '.F90'
    else:
        if not f90_filename.endswith('.f90'):
            f90_filename += '.f90'


    f = open(f90_filename, 'w')

    # Store reaction type: 1 for normal, 2 for thirdbody, 4 for fall-off
    # Store reaction direction: any reaction initially in cti is labeled as forward (reversible or not)

    reac_type = []
    reac_direction = []
    reac_names = []
    for i in range(nr):
        reac_type.append(ctmech.reaction_type(i))
        reac_direction.append('f')
        reac_names.append(ctmech.reaction_equation(i))

    # Store reaction reversibility in is_reversible (1 if reversible, 0 otherwise), size = nr
    # Store reaction index of corresponding reverse reaction in reverse_index, size = nr
    # Extend reac_type by nr_reverse elements (reverse reac has same type as forward)
    # Extend reac_direction by nr_reverse elements (those will all contain "b")

    is_reversible = []
    reverse_index = []
    nr_reverse = 0
    for i in range(nr):
        if ctmech.is_reversible(i):
            is_reversible.append(1)
            reac_type.append(ctmech.reaction_type(i))
            reac_direction.append('b')
            reac_names.append(ctmech.reaction_equation(i))
            reverse_index.append(nr + nr_reverse)
            nr_reverse += 1
        else:
            is_reversible.append(0)
            reverse_index.append(0)

    # Define numerical labels for all reactions. Reverse get same label as corresponding forward

    reac_label = np.zeros(nr + nr_reverse, dtype=int)
    k = 0
    for i in range(nr):
        reac_label[k] = i + 1
        k += 1
    for i in range(nr):
        if is_reversible[i] == 1:
            reac_label[k] = i + 1
            k += 1

    # ThreeBodiesindex: size nTB+nTB_reverse, contains the list of all TB reaction indices (type = 2)
    # nTB: Number of entries in ThreeBodiesindex for forward (or irreversible) reactions
    # nTB_reverse: Number of entries in ThreeBodiesindex for reverse reactions

    # FallOffindex: size nFO+nFO_reverse, contains the list of all FO reaction indices (type = 4)
    # nTB: Number of entries in ThreeBodiesindex for forward (or irreversible) reactions
    # nTB_reverse: Number of entries in ThreeBodiesindex for reverse reactions

    # Plogindex: size nPl+nPl_reverse, contains the list of all Plog reaction indices (type = 5)
    # nPl: Number of entries in Plogindex for forward (or irreversible) reactions
    # nPl_reverse: Number of entries in Plogindex for reverse reactions

    # Initialization
    ThreeBodiesindex = []
    FallOffindex = []
    Plogindex = []
    nTB = 0
    nFO = 0
    nPlog = 0
    nTB_reverse = 0
    nFO_reverse = 0
    nPlog_reverse = 0
    Plog_pressures = []

    # Loop through all reactions, forward and backward
    for i in range(nr + nr_reverse):

        # Simple Three Body reactions
        if reac_type[i] == 2:
            ThreeBodiesindex.append(i)
            if i < nr:
                nTB += 1
                if is_reversible[i]:
                    nTB_reverse += 1

        # Fall Off reactions (e.g. TROE)
        elif reac_type[i] == 4:
            FallOffindex.append(i)
            if i < nr:
                nFO += 1
                if is_reversible[i]:
                    nFO_reverse += 1

        # Plog reactions
        elif reac_type[i] == 5:
            Plogindex.append(i)
            Plog_pressures.append(len(reactions[reac_label[i] - 1].rates))
            if i < nr:
                nPlog += len(reactions[i].rates)
                if is_reversible[i]:
                    nPlog_reverse += len(reactions[i].rates)

    if Plogindex:
        Plog_flag = True
    else:
        Plog_flag = False

    # List species that appear with specified efficiency in third body expressions
    # No need to look at reverse reactions, same thirdbodies as forward
    # Will be tested against QSS species list later

    species_in_thirdbody_expressions = []
    TB_reaction_set = [ireac for ireac in ThreeBodiesindex + FallOffindex if ireac < nr]
    for ireac in TB_reaction_set:
        spec_list = ctmech.reaction(ireac).efficiencies
        for spec in spec_list:
            if spec not in species_in_thirdbody_expressions:
                species_in_thirdbody_expressions.append(spec)

    # Connectivity matrices relating species and reactions in mechanism

    net = networks.Network(mechanism)
    nur = net.nur
    nup = net.nup
    nu = nup - nur

    # Skip if no qss
    if nqss > 0:
        # isqss: array of size ns, True if corresponding spec in species is QSS
        isqss = np.isin(species_names, species_qss_names)
        isqss_og = np.isin(species_names_original, species_qss_names_original)

        # qss_index_global: array of size nqss, containing index of qss species in global species array
        qss_index_global = np.array(np.argwhere(isqss == True))
        qss_index_global_og = np.array(np.argwhere(isqss_og == True))

        # Sort QSS species array according to order of global species array (coming from cti)
        qss_index_local = {}
        qss_sorted = []
        qss_sorted_og = []
        for i in range(nqss):
            qss_sorted.append(species_names[int(qss_index_global[i])])
            qss_sorted_og.append(species_names_original[int(qss_index_global_og[i])])
            qss_index_local[species_names[int(qss_index_global[i])]] = i
        species_qss_names = qss_sorted
        species_qss_names_original = qss_sorted_og

        # QSS/QSS species connectivity

        # Get connectivity between species and reactions
        deltass_qss = net.deltaSS
        deltasr_qss = net.deltaSR

        # Trim to only keep QSS/QSS connectivity
        j = 0
        for i in range(ns):
            if i not in qss_index_global:
                deltass_qss = np.delete(deltass_qss, j, 0)
                deltass_qss = np.delete(deltass_qss, j, 1)
                deltasr_qss = np.delete(deltasr_qss, j, 0)
            else:
                j += 1

        # non zero entries in QSS/QSS connectivity matrix
        index_coupled_qss_i, index_coupled_qss_j = np.nonzero(deltass_qss)
        index_qss_reactions_i, index_qss_reactions_j = np.nonzero(deltasr_qss)

    else:
        qss_index_global = 0
        qss_index_local = 0
        deltass_qss = []
        index_coupled_qss_i = []
        index_coupled_qss_j = []
        index_qss_reactions_i = []
        index_qss_reactions_j = []


    # =================================================================
    # Transfer arrays
    # =================================================================

    # Dictionary for variables
    constants = {}
    mech_variables = {}
    print_variables = {}
    reactions_variables = {}
    qss_variables = {}

    constants['Rcst'] = Rcst
    constants['ne'] = ne
    constants['ns'] = ns
    constants['nr'] = nr
    constants['nr_reverse'] = nr_reverse
    constants['nTB'] = nTB
    constants['nTB_reverse'] = nTB_reverse
    constants['nFO'] = nFO
    constants['nFO_reverse'] = nFO_reverse
    constants['nPlog'] = nPlog
    constants['nPlog_reverse'] = nPlog_reverse
    constants['nqss'] = nqss

    mech_variables['mech'] = ctmech
    mech_variables['elements'] = elements
    mech_variables['n_atoms_in_spec'] = n_atoms_in_spec
    mech_variables['species'] = species_names
    mech_variables['species_og'] = species_names_original
    mech_variables['reactions'] = reactions
    mech_variables['nur'] = nur
    mech_variables['nup'] = nup
    mech_variables['nu'] = nu
    mech_variables['indr'] = net.indr
    mech_variables['inds'] = net.inds
    mech_variables['timescales'] = mechanism.timescales

    reactions_variables['reac_type'] = reac_type
    reactions_variables['reac_label'] = reac_label
    reactions_variables['reac_names'] = reac_names
    reactions_variables['reac_direction'] = reac_direction
    reactions_variables['is_reversible'] = is_reversible
    reactions_variables['ThreeBodiesindex'] = ThreeBodiesindex
    reactions_variables['FallOffindex'] = FallOffindex
    reactions_variables['Plogindex'] = Plogindex
    reactions_variables['Plog_pressures'] = Plog_pressures
    reactions_variables['reverse_index'] = reverse_index

    qss_variables['species_qss_names'] = species_qss_names
    qss_variables['species_qss_names_og'] = species_qss_names_original

    mech_variables['non_qss_species'] = [spec for spec in species_names if spec not in species_qss_names]
    mech_variables['non_qss_species_og'] = \
        [spec for spec in species_names_original if spec not in species_qss_names_original]
    if species_qss_names:
        qss_variables['qss_index_global'] = qss_index_global
        qss_variables['qss_index_local'] = qss_index_local
        qss_variables['deltass_qss'] = deltass_qss
        qss_variables['index_coupled_qss_i'] = index_coupled_qss_i
        qss_variables['index_coupled_qss_j'] = index_coupled_qss_j
        qss_variables['index_qss_reactions_i'] = index_qss_reactions_i
        qss_variables['index_qss_reactions_j'] = index_qss_reactions_j

    # =================================================================
    # Arrhenius parameters
    # =================================================================

    forward_arrh = np.zeros([3, nr])
    forward_arrh_0 = np.zeros([3, nFO])
    forward_arrh_inf = np.zeros([3, nFO])
    forward_arrh_pdep = []

    reverse_arrh = np.zeros([10, nr])
    reverse_arrh_0 = np.zeros([10, nFO])
    reverse_arrh_inf = np.zeros([10, nFO])
    reverse_arrh_pdep = []

    iFO = 0
    for iforw in range(nr):
        # Initialization of lists, cannot be done with matrices as shape is not fixed
        forward_arrh_pdep.append([])
        reverse_arrh_pdep.append([])

        # Arrhenius parameters - distinguish between fall-off and non-fall-off
        # Pre-exponential factors are converted to mol from kmol, with correction factor
        orders = ctmech.reaction(reac_label[iforw] - 1).orders
        if orders:
            reaction_order = sum([orders[spec] for spec in ctmech.reaction(reac_label[iforw] - 1).reactants])
        else:
            reaction_order = sum(nur[:, iforw])

        A_correction_factor = 10 ** (3 * (reaction_order - 1))
        A_with_M_correction_factor = 10 ** (3 * reaction_order)
        # Activation energies are divided by 1000 to convert to J/mol from J/kmol
        E_correction_factor = 1000

        # Multiplier used for sensitivity analysis
        multip = ctmech.multiplier(iforw)

        # Record arrhenius parameters
        if reac_type[iforw] == 4:

            # Get Cantera Arrhenius rates
            my_arrhenius_inf = ctmech.reaction(iforw).high_rate
            my_arrhenius_0 = ctmech.reaction(iforw).low_rate

            # -------------------------------------------
            # Store Arrhenius parameters for forward, fall-off reaction, converted to SI - INF
            forward_arrh_inf[0, iFO] = multip * my_arrhenius_inf.pre_exponential_factor / A_correction_factor
            forward_arrh_inf[1, iFO] = my_arrhenius_inf.temperature_exponent
            forward_arrh_inf[2, iFO] = my_arrhenius_inf.activation_energy / E_correction_factor

            # Store Arrhenius parameters for forward, fall-off reaction, converted to SI - 0
            forward_arrh_0[0, iFO] = multip * my_arrhenius_0.pre_exponential_factor / A_with_M_correction_factor
            forward_arrh_0[1, iFO] = my_arrhenius_0.temperature_exponent
            forward_arrh_0[2, iFO] = my_arrhenius_0.activation_energy / E_correction_factor

            # -------------------------------------------
            # Store Arrhenius parameters for reverse, fall-off reaction, converted to SI
            if is_reversible[iforw] == 1:
                # INF
                reverse_arrh_inf[:, iFO] = reverse_arrhenius_fit(constants, mech_variables, nu[:, iforw],
                                                                 forward_arrh_inf[:, iFO], fit_order=fit_order,
                                                                 start_temperature=start_temperature,
                                                                 end_temperature=end_temperature,
                                                                 number_of_temperatures=number_of_temperatures,
                                                                 pressure=pressure)
                # 0
                reverse_arrh_0[:, iFO] = reverse_arrhenius_fit(constants, mech_variables, nu[:, iforw],
                                                               forward_arrh_0[:, iFO], fit_order=fit_order,
                                                               start_temperature=start_temperature,
                                                               end_temperature=end_temperature,
                                                               number_of_temperatures=number_of_temperatures,
                                                               pressure=pressure)

            # Increment counter for reverse Fall-off reactions
            iFO += 1

        elif reac_type[iforw] == 5:

            # Get Cantera Arrhenius rates
            my_arrhenius = ctmech.reaction(iforw).rates

            for index, arr_set in enumerate(my_arrhenius):

                mylist = []

                # -------------------------------------------
                # Store Arrhenius parameters for forward, non-fall-off reaction, converted to SI
                # Order of reaction is what it is
                mylist.append(arr_set[1].pre_exponential_factor * multip)  # / A_correction_factor)
                mylist.append(arr_set[1].temperature_exponent)
                mylist.append(arr_set[1].activation_energy / E_correction_factor)
                # Storing pressure
                mylist.append(arr_set[0])

                forward_arrh_pdep[iforw].append(mylist)

                # -------------------------------------------
                # Store Arrhenius parameters for reverse, non-fall-off reaction, converted to SI
                if is_reversible[iforw] == 1:
                    reverse_list = list(reverse_arrhenius_fit(constants, mech_variables, nu[:, iforw],
                                                              mylist[:3], fit_order=fit_order,
                                                              start_temperature=start_temperature,
                                                              end_temperature=end_temperature,
                                                              number_of_temperatures=number_of_temperatures,
                                                              pressure=pressure))
                    reverse_list.append(arr_set[0])

                    reverse_arrh_pdep[iforw].append(reverse_list)
        else:

            # Get Cantera Arrhenius rates
            my_arrhenius = ctmech.reaction(iforw).rate

            # -------------------------------------------
            # Store Arrhenius parameters for forward, non-fall-off reaction, converted to SI
            if reac_type[iforw] == 2:
                # Third-body: account for M as species
                forward_arrh[0, iforw] = multip * my_arrhenius.pre_exponential_factor / A_with_M_correction_factor
            else:
                # Order of reaction is what it is
                forward_arrh[0, iforw] = multip * my_arrhenius.pre_exponential_factor / A_correction_factor

            forward_arrh[1, iforw] = my_arrhenius.temperature_exponent
            forward_arrh[2, iforw] = my_arrhenius.activation_energy / E_correction_factor

            # -------------------------------------------
            # Store Arrhenius parameters for reverse, non-fall-off reaction, converted to SI
            if is_reversible[iforw] == 1:
                reverse_arrh[:, iforw] = reverse_arrhenius_fit(constants, mech_variables, nu[:, iforw],
                                                               forward_arrh[:, iforw], fit_order=fit_order,
                                                               start_temperature=start_temperature,
                                                               end_temperature=end_temperature,
                                                               number_of_temperatures=number_of_temperatures,
                                                               pressure=pressure)

    arrhenius_constants = {}

    arrhenius_constants['forward_arrh'] = forward_arrh
    arrhenius_constants['forward_arrh_0'] = forward_arrh_0
    arrhenius_constants['forward_arrh_inf'] = forward_arrh_inf
    arrhenius_constants['forward_arrh_pdep'] = forward_arrh_pdep

    arrhenius_constants['reverse_arrh'] = reverse_arrh
    arrhenius_constants['reverse_arrh_0'] = reverse_arrh_0
    arrhenius_constants['reverse_arrh_inf'] = reverse_arrh_inf
    arrhenius_constants['reverse_arrh_pdep'] = reverse_arrh_pdep

    if use in ['Cantera', 'AVBP', 'YALES2', 'NTMIX']:
        precision = 'pr'

    # elif use in ['YALES2']:
    #     precision = 'WP'

    else:
        logger.warning('Non valid keyword for use ! Must be AVBP, YALES2 or Cantera (default)')
        logger.warning('The f90 file will be generated in Cantera format')
        precision = 'pr'

    # Suffix indicating floating point precisions in all f90 files being generated
    print_variables['precision'] = precision
    # Pressure dependent arrhenius needed or not
    print_variables['Plog'] = Plog_flag
    # QSS call needed or not
    QSS_flag = True if nqss > 0 else False
    print_variables['QSS'] = QSS_flag
    TB_flag = True if nTB > 0 else False
    print_variables['TB'] = TB_flag
    # Routine name inside fortran file
    print_variables['routine_name'] = routine_name
    # File name
    print_variables['f90_filename'] = f90_filename
    # rrate routine
    print_variables['rrate'] = rrate
    # flag for ending module
    EM_flag = False if use in ['YALES2'] else True
    print_variables['EM_flag'] = EM_flag
    # exponential chemistry
    print_variables['exponential'] = exponential
    # classical chemistry
    print_variables['classical'] = classical
    # implicit exponential chemistry
    print_variables['implicit_exponential'] = implicit_exponential
    # Supplementary info for header
    print_variables['supp_info'] = supplement_dict
    # File handle
    print_variables['f'] = f

    # =================================================================
    # Print statements
    # =================================================================

    pr.print_header(print_variables, use)
    if semi_implicit:
        semi_implicit_bool = True
    else:
        semi_implicit_bool = False
    pr.print_declarations(print_variables, constants, mech_variables, qss_variables, reactions_variables, use,
                          semi_implicit_bool)
    if use == 'AVBP':
        pr.print_reaction_expressions(print_variables, constants, reactions_variables)
    if FallOffindex:
        pr.print_lindemann(print_variables)
    if nTB + nFO > 0:
        pr.print_third_body(print_variables, constants, mech_variables, reactions_variables)
    if Plogindex:
        pr.print_pdep(print_variables, constants, arrhenius_constants, mech_variables, reactions_variables)
    pr.print_rate_coeff(print_variables, constants, arrhenius_constants, mech_variables, reactions_variables)
    pr.print_reaction_rates(print_variables, constants, mech_variables, qss_variables, reactions_variables)
    if exponential:
        pr.print_production_rates_exp(print_variables, constants, mech_variables, qss_variables, reactions_variables, use)
        pr.print_gauss_elimination_intern(print_variables)
    else:
        if use != 'NTMIX':
            pr.print_production_rates(print_variables, constants, mech_variables, qss_variables, reactions_variables)
    if use == 'NTMIX' and classical :
         pr.print_production_rates(print_variables, constants, mech_variables, qss_variables, reactions_variables)
    if implicit_exponential:
        pr.print_production_rates_exp(print_variables, constants, mech_variables, qss_variables, reactions_variables, use)
    # Determine and print expressions for qss species if needed
    if species_qss_names:
        qss.print_qss_concentrations(print_variables, constants, mech_variables, reactions_variables, qss_variables,
                                     force=force, maintain_given_order=maintain_given_order)

    if use == 'Cantera':
        pr.print_mass_to_concentration(print_variables,use)
        pr.print_cantera_specific(print_variables, dt_exp_user)

    elif use == 'AVBP':
        pr.print_avbp_specific(print_variables, mech_variables, reactions_variables, semi_implicit)

    elif use == 'YALES2':
        pr.print_mass_to_concentration(print_variables,use)
        pr.print_yales2_specific(print_variables)

    elif use == 'NTMIX':
        pr.print_mass_to_concentration(print_variables,use)
        pr.print_ntmix_specific(print_variables)

    f.close()

    if verbose:
        logger.info(f"Fortran file successfully written as {f90_filename}")


def compile_fortran(f90_filename, force=False, fortran_format='f90',
                    mechanism_lib=None, output=False, remote_install=True):
    """Creates dynamic library to compute chemistry source terms for Cantera.
    Compile a f90 file (potentially generated through create_mech_f90) containing
    all routines necessary to return wdot to Cantera during runs. The output is a
    dynamic library, used in Cantera via the custom kinetics keyword.

    Parameters
    ----------
    f90_filename : str
        name of the f90 file to compile
    force : bool
        if True, compiles fortran even if it hasn't been changed (Default value = False)
    fortran_format : str
        string specifying if it's for example in f77 format (Default value = 'f90')
    mechanism_lib : str
        path to the directory where dynamic libraries are to be compiled (Default value = None)
    output : bool
        if True displays the output of the make (Default value = False)
    remote_install : bool
        if False the compilation will be done automatically in the sources (Default value = True)


    """

    ## Test to know if subroutine has the right name to do the compiling
    line_custom = None
    f = open(f90_filename, 'r')
    for lines in f.readlines():
        if 'subroutine' in lines.lower() and '(p,t,y,wdot)' in lines.replace(" ", "").lower():
            line_custom = lines
            break
    f.close()

    if line_custom == None:
        print("The main subroutine of the f90 file should be written : 'subroutine customkinetics (P, T, y, wdot)' to work")
        sys.exit()
    else:
        if 'customkinetics' not in line_custom:
            print("The main subroutine of the f90 file should be written : 'subroutine customkinetics (P, T, y, wdot)' to work")
            sys.exit()

    if not remote_install:
        # Finds the path to the cantera module and then copies the f90 file in the correct directory
        path_to_init = ct.__file__

        terminator = path_to_init.rfind('/lib')
        path = (path_to_init[:terminator])

        lib_name = '/lib'

        path_to_dir = path + lib_name

        if not os.path.isdir(path_to_dir) and os.path.isdir(path_to_dir + '64'):
            lib_name = '/lib64'
            path_to_dir = path + lib_name

        elif not os.path.isdir(path_to_dir) and not os.path.isdir(path_to_dir + '64'):
            logger.error('There is a problem with your installation, ' + path_to_dir + ' does not exist')
            sys.exit()

    elif 'CUSTOM_LIB' not in os.environ:
        if mechanism_lib:
            path_to_dir = mechanism_lib
        else:
            path_to_file = os.path.realpath(__file__)
            terminator = path_to_file.rfind('/')
            path = (path_to_file[:terminator])
            lib_name = '/mech_lib'

            path_to_dir = path + lib_name

        if not os.path.isdir(path_to_dir):
            os.mkdir(path_to_dir)
            head, tail = display.head_tail_charac(text='WARNING', size='small')
            logger.info(head)
            logger.info('The directory ' + path_to_dir + ' has been created')
            logger.info('The necessary files for the compilation and customised kinetics run will be stored there')
            logger.info("If this location does not suit you, please specify a path as mechanism_lib argument")
            logger.info(tail)

    else:
        path_to_dir = os.environ['CUSTOM_LIB']

    # Correct path definitions

    path_to_ext = path_to_dir + '/customkinetics.f90'
    path_to_makefile = path_to_dir + '/Makefile'

    make_call = 'make -C ' + path_to_dir
    if not output:
        make_call += ' &> compilation_output.txt'

    # Copying f90 file to the right place

    if not f90_filename.endswith('.f90'):
        f90_filename = f90_filename + '.f90'

    if not os.path.isfile(f90_filename):
        if os.path.isfile(tools.database + '/mechs/' + f90_filename):
            f90_filename = tools.database + '/mechs/' + f90_filename

    if not os.path.isfile(path_to_ext) or not filecmp.cmp(f90_filename, path_to_ext) or force:

        if not os.path.isdir(path_to_dir):
            logger.info(path_to_dir + 'directory created')
            database.create_dir(path_to_dir)

        # logger.info(path_to_ext)
        shutil.copy(f90_filename, path_to_ext)

        # Checks if the Makefile exists, if not creates it

        if not os.path.isfile(path_to_makefile):

            f = open(path_to_makefile, 'w')
            text = """\
customkinetics.so: customkinetics.f90
\tgfortran -ffixed-line-length-0 -c customkinetics.f90 -g -fPIC -o customkinetics.o
\tgfortran -shared -o customkinetics.so customkinetics.o
    """
            if fortran_format == 'f77':
                text = """\
customkinetics.so: customkinetics.f90
\tgfortran -ffixed-line-length-0 -ffixed-form -c customkinetics.f90 -g -fPIC -o customkinetics.o
\tgfortran -shared -o customkinetics.so customkinetics.o
    """
            f.write(text)
            f.close()
            #head, tail = display.head_tail_charac(text='WARNING', size='small')
            #logger.info(head)
            logger.info('The Makefile was missing or wrongly named')
            logger.info('The correct Makefile has been created (' + path_to_makefile + ')')
            #logger.info(tail)

    # In remote installation, sets the correct LD_LIBRARY_PATH for the dynamic library
    if 'CUSTOM_LIB' not in os.environ:
        custom_env_set = 'export CUSTOM_LIB=' + path_to_dir
        environment_set = 'export LD_LIBRARY_PATH=$CUSTOM_LIB:$LD_LIBRARY_PATH'
        if sys.platform == 'darwin':
            environment_set += '\nexport DYLD_LIBRARY_PATH=$CUSTOM_LIB:$DYLD_LIBRARY_PATH'

        if mechanism_lib:
            dir_status = 'you chose'
            text = ''
        else:
            dir_status = 'automatically created'

            text = ("""\
This part is not automated as it would require sneaky modifications of your bashrc and that is privacy violation !
ARCANE is not that kind of code !

Please copy those commands in your bashrc (or execute them in your shell):
{0}
{1}

It will add the directory {2} to the dynamic library path blab bla bla, informatics stuff ...

I advise you stick to only one directory for the compilation for 2 reasons:
- 1: it will be boring to add the path again.
- 2: it's absolutely useless to have several.

That's it ! Enjoy using ARCANE.""")

        head, tail = display.head_tail_charac(text='WARNING', size='medium', style='#')

        logger.critical(head)
        display.print_greetings_and_farewell(greeting=True)
        logger.critical(text.format(custom_env_set, environment_set, dir_status))
        display.print_greetings_and_farewell(farewell=True)
        logger.critical(tail)
        sys.exit()

    process = subprocess.call(make_call, shell=True)

    if process != 0:
        logger.error('The compilation of the f90 file failed ! There is something wrong ...')
        if not output:
            logger.error('Check compilation_output.txt for intel')
        exit()
    else:

        logger.info('Compilation of the f90 successful')
        subprocess.call('rm -f compilation_output.txt', shell=True)
        time.sleep(0.05)  # Necessary for python to understand there was a change


def cantera_qss_solution(ctmech, species_qss_names, qss_cti_file=None, transport=None, output_solution=True):
    """Creates a Cantera solution object from the provided cti, with the specified list of QSS species.

    Parameters
    ----------
    ctmech :
        `Cantera.Solution` object
    species_qss_names : list
        list of qss species names
    qss_cti_file : str
        name of cti file to print with the proper list of qss species (Default value = None)
    transport : str
        transport model (Default value = None)
    output_solution : bool
        if False, the Solution object will not be loaded (Default value = True)

    Returns
    -------
    phase:
        Cantera solution object with the right species and reactions and with custom kinetics

    """
    ct.suppress_thermo_warnings()

    # Species objects list
    species_objects = ctmech.species()

    # Species list
    species_names = ctmech.species_names

    # Species objects list sorted according to the species list
    species_objects_names = []
    for spec in species_objects:
        species_objects_names.append(spec.name)

    # Corrects inconstistency between the list of species and the list of species objects
    if species_objects_names != species_names:
        species_objects_sorted = []
        for spec in species_names:
            species_objects_sorted.append(species_objects[species_objects_names.index(spec)])
    else:
        species_objects_sorted = species_objects

    # Species non QSS
    species_non_qss = [s for s in species_objects_sorted if s.name not in species_qss_names]
    species_non_qss_names = [s.name for s in species_non_qss]

    transport_model = ctmech.transport_model

    fuel = ''
    H_atoms = 0

    if 'H' not in ctmech.element_names:
        global_reaction = "reaction('" + ctmech.species_names[0] + " => " \
                          + ctmech.species_names[0] \
                          + "', [0.0, 0.0, 0.0])"
    else:
        for spec in species_non_qss_names:
            if ctmech.n_atoms(spec, 'H') > H_atoms:
                fuel = spec
                H_atoms = ctmech.n_atoms(spec, 'H')

        if 'O2' in species_non_qss_names \
                and 'CO2' in species_non_qss_names \
                and 'H2O' in species_non_qss_names :
            nC = ctmech.n_atoms(fuel, 'C') if 'C' in ctmech.species(fuel).composition else 0
            nH = ctmech.n_atoms(fuel, 'H') if 'H' in ctmech.species(fuel).composition else 0
            nO = ctmech.n_atoms(fuel, 'O') if 'O' in ctmech.species(fuel).composition else 0

            stoich = nC + 0.25 * nH - 0.5 * nO

            if nC != 0 and 'CO2' in species_non_qss_names:
                CO2string = str(nC) + " CO2 + "
                global_reaction = "reaction('" + fuel + " + " + str(stoich) + " O2 => " \
                                  + CO2string + str(nH / 2) + " H2O " \
                                  + "', [0.0, 0.0, 0.0])"
            elif nC != 0 and 'CO2' not in species_non_qss_names:
                global_reaction = "reaction('" + fuel + " => " \
                                  + fuel \
                                  + "', [0.0, 0.0, 0.0])"
            else:
                CO2string = ""
                global_reaction = "reaction('" + fuel + " + " + str(stoich) + " O2 => " \
                                  + CO2string + str(nH / 2) + " H2O " \
                                  + "', [0.0, 0.0, 0.0])"

        else:
            global_reaction = "reaction('" + fuel + " => " \
                              + fuel \
                              + "', [0.0, 0.0, 0.0])"

    # Generating the cti file (if qss_cti_file= None, a dummy one will be generated)
    if not qss_cti_file:
        qss_cti_file = "dummy.cti"

    ct2cti.write(ctmech, qss_cti_file, kinetics='custom', transport=transport,
                 species_qss=species_qss_names, dummy=global_reaction)

    if output_solution:
        phase = ct.Solution(qss_cti_file)

    if qss_cti_file == "dummy.cti":
        subprocess.call('rm -rf ' + qss_cti_file, shell=True)

    if output_solution:
        return phase


def cantera_reduced_solution(ctmech, species_to_discard=[], reactions_to_discard=[],
                             diluent=[], auto_remove=False, reduced_cti_file=None,
                             transport=None):
    """Removes species or reactions from a Cantera Solution object.

    Parameters
    ----------
    ctmech :
        `Cantera.Solution` object
    species_to_discard : list
        list of species to discard (Default value = [])
    reactions_to_discard : list
        list of reactions indices to discard (Default value = [])
    diluent : list
        species that even inactive will be kept (Default value = [])
    auto_remove : bool
        if True, removes species that are no longer used in a reaction (Default value = False)
    reduced_cti_file : str
        name of the output .cti (Default value = None)
    transport : str
        transport model (Default value = None)

    Returns
    -------
    modified Cantera Solution object

    """
    ct.suppress_thermo_warnings()
    # Species objects list
    species_objects = ctmech.species()

    # Species list
    species_names = ctmech.species_names

    # Reactions objects list
    reactions_objects = ctmech.reactions()

    reactions_to_discard_equations = [reac.equation for reac in reactions_to_discard]

    # Initialization
    reactions_to_include = []
    # Species objects list sorted according to the species list
    species_objects_names = []
    for spec in species_objects:
        species_objects_names.append(spec.name)

    # Corrects inconsistency between the list of species and the list of species objects

    if species_objects_names != species_names:
        species_objects_sorted = []
        for spec in species_names:
            species_objects_sorted.append(species_objects[species_objects_names.index(spec)])
    else:
        species_objects_sorted = species_objects

    # Species to include in new mechanism
    species_to_include = [s for s in species_objects_sorted if s.name not in species_to_discard]
    species_to_include_names = [s.name for s in species_to_include]

    species_in_a_reaction = set()

    # Reactions to include in new mechanism
    for id_reac, reac in enumerate(reactions_objects):

        if not all(reactant in species_to_include_names for reactant in reac.reactants):
            continue

        if not all(product in species_to_include_names for product in reac.products):
            continue

        # Species presents in a reaction
        species_in_a_reaction.update(reac.reactants)
        species_in_a_reaction.update(reac.products)

        if reac.reaction_type in [2, 4]:
            new_efficiencies = reac.efficiencies
            species_in_a_reaction.update(new_efficiencies.keys())
            for spec in species_to_discard:
                if spec in reac.efficiencies:
                    new_efficiencies.pop(spec)
            reac.efficiencies = new_efficiencies

        if reac.equation in reactions_to_discard_equations:
            if reac.duplicate:
                i = reactions_to_discard_equations.index(reac.equation)
                reac_discard = reactions_to_discard[i]
                if reac.reaction_type == 5:
                    if reac.rates == reac_discard.rates:
                            continue
                elif reac.reaction_type == 4:
                    if reac.high_rate == reac_discard.high_rate and reac.low_rate == reac_discard.low_rate:
                        continue
                else:
                    if reac.rate == reac_discard.rate:
                        continue
            else:
                continue

        reactions_to_include.append(reac)

    species_in_a_reaction.update(diluent)
    # If a species is no more in a reaction, exclude it
    automatically_removed = list(set(species_to_include_names) - set(species_in_a_reaction))
    if automatically_removed and species_to_discard and auto_remove:
        logger.warning('WARNING: Species ' + ', '.join(
            automatically_removed) + ' automatically removed since it was no longer present in a reaction')
        species_to_include = [s for s in species_to_include if s.name in species_in_a_reaction]

    # Generating the cti file (if qss_cti_file= None, a dummy one will be generated)
    if not reduced_cti_file:
        reduced_cti_file = "dummy.cti"

    phase = ct.Solution(thermo='IdealGas',
                        transport=transport,
                        kinetics='GasKinetics',
                        species=species_to_include,
                        reactions=reactions_to_include)

    ct2cti.write(phase, reduced_cti_file, transport=transport)

    if reduced_cti_file == "dummy.cti":
        subprocess.call('rm -rf ' + reduced_cti_file, shell=True)

    return phase


def reverse_arrhenius_fit(constants, mech_variables, nu, forward_arrh, fit_order=0,
                          start_temperature=300, end_temperature=3000, number_of_temperatures=6000,
                          pressure=101325):
    r"""Computes the Arrhenius coefficients of the backward reactions by computing the reverse rate constants
    with the equilibrium constants on a set of Temperature. The coefficients are then fitted with the
    `numpy.opt.curve_fit` function

    Parameters
    ----------
    constants : dict
        dictionary of constants needed for computation
    mech_variables : dict
        dictionary of mechanism related parameters needed for computation
    nu : list
        net stoichiometric coefficients :math:`\nu` for reaction (size :math:`N_s`)
    forward_arrh : list
        Arrhenius coefficients of the reversible forward reaction, according to [:math:`A`, :math:`\beta`, :math:`E_a`]
    fit_order : int
        order of the fit used for the modified Arrhenius function (Default value = 0)
    start_temperature : float
        start temperature for the fit (Default value = 300)
    end_temperature : float
        end temperature for the fit (Default value = 3000)
    number_of_temperatures : int
        number of temperature used for the fit (Default value = 6000)
    pressure : float
        pressure used in the equilibrium constant computation (Default value = 101325)


    Can be derived from [1]_, with :math:`\ln\left(K_j(T)\right) = - \sum_{i=N_s} \nu_j \mu_j^0`
    and :math:`\mu_i^0 = \mu_i - RT \ln\left(\frac{P_i}{P_0}\right)`, leading to:

    .. math::
            \ln\left(K_j(T)\right) = - \frac{1}{RT}  \sum_{i=N_s} \nu_{i,j} \left(\mu_i - RT\ln\left(\frac{P_i}{P_0}\right)\right) \\\\
            = - \frac{1}{RT}  \sum_{i=N_s} \nu_{i,j} \mu_i - \sum_{i=1,N_s} \nu_i \ln\left(\frac{P_i}{P_0}\right)

    with, :math:`P_i = \frac{\rho Y_i}{W_i}RT \approx RT\quad` if   :math:`Y_i = 1`, :math:`W_i = \rho`, finally:

    .. math::
            \left.\ln\left(K_j(T)\right)\right|_P = \sum\nu_j \ln\left(\frac{RT}{P_0}\right) - \sum\frac{\nu_j \mu_j}{RT}


    .. [1] N. Peters "*Fifteen Lectures on Laminar and Turbulent Combustion*" (1992) p. 14 eq.1.71

    Returns
    -------
    reverse_Arrh : np.array
        Array of Arrhenius coefficients (:math:`A`, :math:`\beta` and :math:`E_a`) of the corresponding backward reaction

    """

    Rcst = constants['Rcst']

    ctmech = mech_variables['mech']

    # Temperature array, equidistant in 1/T
    temperatures = np.zeros([number_of_temperatures])
    dx = (1.0 / end_temperature - 1.0 / start_temperature) / (number_of_temperatures - 1.0)
    for i in range(number_of_temperatures):
        temperatures[i] = 1.0 / (1.0 / start_temperature + dx * i)

    RT_P = (Rcst * temperatures) / pressure

    # Non zero indices of nu matrix
    indj = np.nonzero(nu)

    # Equilibrium constant computation for each temperature

    ln_K = np.zeros([number_of_temperatures])
    for iT in range(number_of_temperatures):
        ctmech.TP = temperatures[iT], pressure
        sum_nu = 0
        sum_nu_mu = 0
        for j in indj[0]:
            sum_nu += nu[j]
            sum_nu_mu += nu[j] * ctmech.standard_gibbs_RT[j]

        ln_K[iT] = - sum_nu * np.log(RT_P[iT]) - sum_nu_mu

    # Forward rate constant computation using Arrhenius coefficients
    # Backward rate constant computation using forward and equilibrium
    # Fitting data with numpy function and storing into matrix of Arrhenius coefficients

    A = forward_arrh[0]
    beta = forward_arrh[1]
    Ea = forward_arrh[2]

    if A < 0:
        negative_A = True
        A = - A
    else:
        negative_A = False

    ln_kf = np.log(A) + beta * np.log(temperatures) + (-Ea / (Rcst * temperatures))
    ln_kb = ln_kf - ln_K

    reverse_Arrh, pcov = opt.curve_fit(getattr(sys.modules[__name__], f'calculate_lnk_{fit_order}'),
                                       temperatures, ln_kb)

    debug = False  # True
    if debug:
        import matplotlib.pyplot as plt
        lnb = getattr(sys.modules[__name__], f'calculate_lnk_{fit_order}')(temperatures, *reverse_Arrh)
        plt.plot(np.exp(lnb), 'r-')
        plt.plot(np.exp(ln_kb), 'b--')
        plt.show()

    reverse_Arrh[0] = np.exp(reverse_Arrh[0])

    if negative_A:
        reverse_Arrh[0] = - reverse_Arrh[0]

    reverse_Arrh[2] = reverse_Arrh[2] * Rcst

    reverse_Arrh = list(reverse_Arrh)
    reverse_Arrh += [0] * (7 - fit_order)

    reverse_Arrh = np.array(reverse_Arrh)

    return reverse_Arrh


def calculate_lnk_0(T, ln_A, beta, Ea_R):
    r"""Formula of :math:`ln(k)` using Arrhenius expression, use SI units

    Parameters
    ----------
    T :
        temperature :math:`T`
    ln_A :
        natural log of pre-exponential factor :math:`\ln(A)`
    beta :
        temperature exponent :math:`\beta`
    Ea_R :
        activation energy divided by ideal gas constant, :math:`\frac{E_a}{R}`

    Returns
    -------
    computed `lnk_0` according to :math:`\ln(A) + \beta\ln(T) + -\frac{E_a}{RT}`

    """
    return ln_A + (beta * np.log(T)) + (-Ea_R / T)


def calculate_lnk_1(T, ln_A, beta, Ea_R, alpha1):
    r"""Formula of :math:`ln(k)` using Arrhenius expression, use SI units

    Parameters
    ----------
    T :
        temperature :math:`T`
    ln_A :
        natural log of pre-exponential factor :math:`\ln(A)`
    beta :
        temperature exponent :math:`\beta`
    Ea_R :
        activation energy divided by ideal gas constant, :math:`\frac{E_a}{R}`
    alpha1 :
        fit function coefficient :math:`\alpha_1`

    Returns
    -------
    computed `lnk_1` according to :math:`\ln(A) + \beta\ln(T) + -\frac{E_a}{RT} + \alpha_1 T`

    """
    return ln_A + (beta * np.log(T)) + (-Ea_R / T) + (alpha1 * T)


def calculate_lnk_2(T, ln_A, beta, Ea_R, alpha1, alpha2):
    r"""Formula of :math:`ln(k)` using Arrhenius expression, use SI units

    Parameters
    ----------
    T :
        temperature :math:`T`
    ln_A :
        natural log of pre-exponential factor :math:`\ln(A)`
    beta :
        temperature exponent :math:`\beta`
    Ea_R :
        activation energy divided by ideal gas constant :math:`\frac{E_a}{R}`
    alpha1 :
        fit function coefficient :math:`\alpha_1`
    alpha2 :
        fit function coefficient :math:`\alpha_2`

    Returns
    -------
    computed `lnk_2` according to :math:`\ln(A) + \beta\ln(T) + -\frac{E_a}{RT} + \alpha_1 T+ \alpha_2 T^2`

    """
    return ln_A + (beta * np.log(T)) + (-Ea_R / T) + (alpha1 * T) + (alpha2 * T ** 2)


def calculate_lnk_3(T, ln_A, beta, Ea_R, alpha1, alpha2, alpha3):
    r"""Formula of :math:`ln(k)` using Arrhenius expression, use SI units

    Parameters
    ----------
    T :
        temperature :math:`T`
    ln_A :
        natural log of pre-exponential factor :math:`\ln(A)`
    beta :
        temperature exponent :math:`\beta`
    Ea_R :
        activation energy divided by ideal gas constant :math:`\frac{E_a}{R}`
    alpha1 :
        fit function coefficient :math:`\alpha_1`
    alpha2 :
        fit function coefficient :math:`\alpha_2`
    alpha3 :
        fit function coefficient :math:`\alpha_3`

    Returns
    -------
    computed `lnk_3` according to :math:`\ln(A) + \beta\ln(T) + -\frac{E_a}{RT} + \alpha_1 T+ \alpha_2 T^2 + \alpha_3 T^3`

    """
    return ln_A + (beta * np.log(T)) + (-Ea_R / T) + (alpha1 * T) + (alpha2 * T ** 2) \
           + (alpha3 * T ** 3)


def calculate_lnk_4(T, ln_A, beta, Ea_R, alpha1, alpha2, alpha3, alpha4):
    r"""Formula of :math:`ln(k)` using Arrhenius expression, use SI units

    Parameters
    ----------
    T :
        temperature :math:`T`
    ln_A :
        natural log of pre-exponential factor :math:`\ln(A)`
    beta :
        temperature exponent :math:`\beta`
    Ea_R :
        activation energy divided by ideal gas constant :math:`\frac{E_a}{R}`
    alpha1 :
        fit function coefficient :math:`\alpha_1`
    alpha2 :
        fit function coefficient :math:`\alpha_2`
    alpha3 :
        fit function coefficient :math:`\alpha_3`
    alpha4 :
        fit function coefficient :math:`\alpha_4`

    Returns
    -------
    computed `lnk_4` according to :math:`\ln(A) + \beta\ln(T) + -\frac{E_a}{RT} + \alpha_1 T+ \alpha_2 T^2 + \alpha_3 T^3 + \alpha_4 T^4`

    """
    return ln_A + (beta * np.log(T)) + (-Ea_R / T) + (alpha1 * T) + (alpha2 * T ** 2) \
           + (alpha3 * T ** 3) + (alpha4 * T ** 4)


def set_qss_coupling_variables(mechanism, species_qss_names):
    """Computes necessary data for qss coupling check

    Parameters
    ----------
    mechanism :
        `ARCANE.mechanisms.Mechanism` object that contains all mechanism-specific info
    species_qss_names : list
        list of the QSS species, can be empty.

    Returns
    -------
    mech_variables : dict
        dictionary containing mechanism information
    reactions_variables : dict
        dictionary containing reaction information
    qss_variables : dict
        dictionary containing QSS information

    """

    ctmech = mechanism.ctmech
    species_names = ctmech.species_names
    reactions = ctmech.reaction_equations()
    ns = len(species_names)
    nr = len(reactions)
    nqss = len(species_qss_names)

    # isqss: array of size ns, True if corresponding spec in species is QSS
    isqss = np.isin(species_names, species_qss_names)

    # qss_index_global: array of size nqss, containing index of qss species in global species array
    qss_index_global = np.array(np.argwhere(isqss == True))

    # Sort QSS species array according to order of global species array (coming from cti)
    qss_index_local = {}
    qss_sorted = []
    for i in range(nqss):
        qss_sorted.append(species_names[int(qss_index_global[i])])
        qss_index_local[species_names[int(qss_index_global[i])]] = i
    species_qss_names = qss_sorted

    net = networks.Network(mechanism)
    nur = net.nur
    nup = net.nup

    is_reversible = []
    for i in range(nr):
        if ctmech.is_reversible(i):
            is_reversible.append(1)
        else:
            is_reversible.append(0)

    # isqss: array of size ns, True if corresponding spec in species is QSS
    isqss = np.isin(species_names, species_qss_names)

    # qss_index_global: array of size nqss, containing index of qss species in global species array
    qss_index_global = np.array(np.argwhere(isqss == True))

    # Get connectivity between species and reactions
    deltass_qss = net.deltaSS
    deltasr_qss = net.deltaSR

    # Trim to only keep QSS/QSS connectivity
    j = 0
    for i in range(ns):
        if i not in qss_index_global:
            deltass_qss = np.delete(deltass_qss, j, 0)
            deltass_qss = np.delete(deltass_qss, j, 1)
            deltasr_qss = np.delete(deltasr_qss, j, 0)
        else:
            j += 1

    # non zero entries in QSS/QSS connectivity matrix
    index_coupled_qss_i, index_coupled_qss_j = np.nonzero(deltass_qss)
    index_qss_reactions_i, index_qss_reactions_j = np.nonzero(deltasr_qss)

    mech_variables = {}
    reactions_variables = {}
    qss_variables = {}

    mech_variables['nur'] = nur
    mech_variables['nup'] = nup

    reactions_variables['is_reversible'] = is_reversible

    qss_variables['species_qss_names'] = species_qss_names
    qss_variables['qss_index_global'] = qss_index_global
    qss_variables['index_coupled_qss_i'] = index_coupled_qss_i
    qss_variables['index_coupled_qss_j'] = index_coupled_qss_j
    qss_variables['index_qss_reactions_i'] = index_qss_reactions_i
    qss_variables['index_qss_reactions_j'] = index_qss_reactions_j

    return mech_variables, reactions_variables, qss_variables


def print_fortran_debug(f90_filename):
    """Prints a f90 file allowing to see if the compilation is good

    Parameters
    ----------
    f90_filename : str
        name of the output .f90 file

    """

    # Establish handle to fortran source file
    f = open(f90_filename + '.f90', 'w')
    text = """\
module mod_customkinetics
  implicit none

end module mod_customkinetics

subroutine customkinetics(P, T, rho, y, wdot)

  use mod_customkinetics
  implicit none

    integer :: nspec = 42
    integer :: nqss = 7
    real, dimension(42) :: y, wdot
    real :: T,P,rho

    write(*,*) "Hello there"
    write(*,*) 'nspec = ', nspec, 'nqss = ', nqss

  return
end subroutine customkinetics
"""

    f.write(text)


def problematic_qss(mechanism, species_qss_list, remove_from='list'):
    """Checks if a QSS species implies quadratic coupling
    or if it's part of the weighted third body species.
    Given the parameter it will discard the species from the list
    or delete the corresponding reactions

    Parameters
    ----------
    mechanism :
        object of class `ARCANE.mechanisms.Mechanism` object
    species_qss_list : list
        list of QSS species to check
    remove_from : str
        if 'list' discard it from the list, while if 'mechanism' delete from the reactions (Default value = 'list')

    """

    # Mechanism variables
    ctmech = mechanism.ctmech
    species_names = ctmech.species_names
    reactions = ctmech.reactions()
    ns = len(species_names)
    nr = len(reactions)

    species_qss_list = tools.sort_species_by_solution(species_qss_list, ctmech)

    # Network variables
    net = networks.Network(mechanism)
    nur = net.nur
    nup = net.nup

    # Initialization of outputs
    new_ctmech = ctmech
    species_qss_names_proposition = species_qss_list

    # Check that all species in QSS list are indeed valid species in mechanism
    qss_not_found = [spec for spec in species_qss_list if spec not in species_names]
    if qss_not_found:
        logger.error('Species ' + str(', ').join(qss_not_found) + ' is/are not in mechanism object')  # + cti_filename)
        sys.exit()

    # Store reaction type: 1 for normal, 2 for thirdbody, 4 for fall-off
    reac_type = []
    for i in range(nr):
        reac_type.append(ctmech.reaction_type(i))
    is_reversible = []
    nr_reverse = 0
    for i in range(nr):
        if ctmech.is_reversible(i):
            is_reversible.append(1)
            reac_type.append(ctmech.reaction_type(i))
            nr_reverse += 1
        else:
            is_reversible.append(0)

    # Initialization
    ThreeBodiesindex = []
    FallOffindex = []
    nTB = 0
    nFO = 0
    nTB_reverse = 0
    nFO_reverse = 0

    # Loop through all reactions, forward and backward
    for i in range(nr + nr_reverse):
        # Simple Three Body reactions
        if reac_type[i] == 2:
            ThreeBodiesindex.append(i)
            if i < nr:
                nTB += 1
                if is_reversible[i]:
                    nTB_reverse += 1

        # Fall Off reactions (e.g. TROE)
        elif reac_type[i] == 4:
            FallOffindex.append(i)
            if i < nr:
                nFO += 1
                if is_reversible[i]:
                    nFO_reverse += 1

    # List species that appear with specified efficiency in third body expressions
    # No need to look at reverse reactions, same thirdbodies as forward
    # Will be tested against QSS species list later

    species_in_thirdbody_expressions = []
    TB_reaction_set = [ireac for ireac in ThreeBodiesindex + FallOffindex if ireac < nr]
    for ireac in TB_reaction_set:
        spec_list = ctmech.reaction(ireac).efficiencies
        for spec in spec_list:
            if spec not in species_in_thirdbody_expressions:
                species_in_thirdbody_expressions.append(spec)

    # Check if QSS species are present as specific third bodies*
    spec_in_thirdbody_expressions = list(set(species_qss_list).intersection(species_in_thirdbody_expressions))

    # isqss: array of size ns, True if corresponding spec in species is QSS
    isqss = np.isin(species_names, species_qss_list)

    # qss_index_global: array of size nqss, containing index of qss species in global species array
    qss_index_global = np.array(np.argwhere(isqss == True))

    # Get connectivity between species and reactions
    deltass_qss = net.deltaSS
    deltasr_qss = net.deltaSR

    # Trim to only keep QSS/QSS connectivity
    j = 0
    for i in range(ns):
        if i not in qss_index_global:
            deltass_qss = np.delete(deltass_qss, j, 0)
            deltass_qss = np.delete(deltass_qss, j, 1)
            deltasr_qss = np.delete(deltasr_qss, j, 0)
        else:
            j += 1

    # non zero entries in QSS/QSS connectivity matrix
    index_coupled_qss_i, index_coupled_qss_j = np.nonzero(deltass_qss)
    index_qss_reactions_i, index_qss_reactions_j = np.nonzero(deltasr_qss)

    species_qss_names_proposition = species_qss_list.copy()
    problematic_species = {}
    problematic_species_self = {}
    problematic_reactions = []

    for i in range(len(species_qss_list)):
        # Species coupled with an other one
        for j in index_coupled_qss_j[index_coupled_qss_i == i]:
            if j != i:
                for k in index_qss_reactions_j[index_qss_reactions_i == i]:
                    if (nur[qss_index_global[i], k] >= 1 and nur[qss_index_global[j], k] >= 1) or \
                            ((nup[qss_index_global[i], k] >= 1 and nup[qss_index_global[j], k] >= 1) and
                             is_reversible[k] == 1):
                        if k not in problematic_reactions:
                            problematic_reactions.append(k)
                            logger.debug(
                                '!!! WARNING !!! No quadratic coupling allowed (species ' + species_qss_list[i]
                                + ' and ' + 'species ' + species_qss_list[j] + ' in reaction ' + str(k + 1) + ')')
                        # Dictionnary on how many times a species is involved in a conflict
                        if not species_qss_names_proposition[i] in problematic_species:
                            problematic_species[species_qss_list[i]] = 0
                        else:
                            problematic_species[species_qss_list[i]] += 1

                        if not species_qss_names_proposition[j] in problematic_species:
                            problematic_species[species_qss_list[j]] = 0
                        else:
                            problematic_species[species_qss_list[j]] += 1

        # Species coupled with themselves
        for k in index_qss_reactions_j[index_qss_reactions_i == i]:
            if nur[qss_index_global[i], k] > 1 or (nup[qss_index_global[i], k] > 1 and is_reversible[k] == 1):
                if k not in problematic_reactions:
                    problematic_reactions.append(k)
                    logger.debug('!!! WARNING !!! No quadratic coupling allowed (species ' + species_qss_list[i]
                                   + ' with coefficient ' + str(
                        int(nur[qss_index_global[i], k])) + ' in reaction ' + str(k + 1) + ')')
                # Dictionary on how many times a species is involved in a conflict
                problematic_species_self[species_qss_list[i]] = -1

    # Species in third body expression
    if spec_in_thirdbody_expressions:
        logger.debug('!!! WARNING !!! Invalid QSS species ' + str(', ').join(
            spec_in_thirdbody_expressions) + ' in thirdbody expressions')

    if remove_from == 'mechanism':

        # If a qss species is present in third body reaction, remove it
        reactions_objects = ctmech.reactions()

        problematic_reactions_objects = [reactions_objects[i] for i in problematic_reactions]

        new_ctmech = cantera_reduced_solution(mechanism.ctmech, reactions_to_discard=problematic_reactions_objects)

        reactions_objects = new_ctmech.reactions()
        reactions_to_include = []

        flag = False
        for reac in reactions_objects:
            if reac.reaction_type in [2, 4]:
                new_efficiencies = reac.efficiencies
                for spec in spec_in_thirdbody_expressions:
                    if spec in reac.efficiencies:
                        new_efficiencies.pop(spec)
                        flag = True
                reac.efficiencies = new_efficiencies

            reactions_to_include.append(reac)

        new_ctmech = ct.Solution(thermo='IdealGas',
                                 transport='Mix',
                                 kinetics='GasKinetics',
                                 species=ctmech.species(),
                                 reactions=reactions_to_include)

        logger.debug("Efficiencies have been updated so they no longer contain QSS species")

        if problematic_reactions:
            problematic_reactions_names = [str(name + 1) for name in problematic_reactions]
            logger.debug('The given subset was not valid due to quadratic coupling'
                         ' and/or QSS species as specific collision partner')
            logger.debug('Implicated reactions and/or efficiencies ' + ', '.join(
                problematic_reactions_names) + ' have been removed')

    else:
        # Proposition of a better QSS set
        while problematic_species:
            species_qss_names_proposition.remove(
                max(problematic_species.keys(), key=(lambda key: problematic_species[key])))
            problematic_species = {}

            for i in range(len(species_qss_list)):
                for j in index_coupled_qss_j[index_coupled_qss_i == i]:
                    if j != i \
                            and species_qss_list[i] in species_qss_names_proposition \
                            and species_qss_list[j] in species_qss_names_proposition:
                        for k in index_qss_reactions_j[index_qss_reactions_i == i]:
                            if (nur[qss_index_global[i], k] >= 1 and nur[qss_index_global[j], k] >= 1) or \
                                    ((nup[qss_index_global[i], k] >= 1 and nup[qss_index_global[j], k] >= 1) and
                                     is_reversible[k] == 1):
                                # Dictionnary on how many times a species is involved in a conflict
                                if species_qss_list[i] not in problematic_species_self:
                                    if not species_qss_list[i] in problematic_species:
                                        problematic_species[species_qss_list[i]] = 0
                                    else:
                                        problematic_species[species_qss_list[i]] += 1

                                if species_qss_list[j] not in problematic_species_self:
                                    if not species_qss_list[j] in problematic_species:
                                        problematic_species[species_qss_list[i]] = 0
                                    else:
                                        problematic_species[species_qss_list[j]] += 1

        for spec in spec_in_thirdbody_expressions:
            if spec in species_qss_names_proposition:
                species_qss_names_proposition.remove(spec)

        for spec in problematic_species_self:
            if spec in species_qss_names_proposition:
                species_qss_names_proposition.remove(spec)

        if species_qss_names_proposition != species_qss_list:
            logger.debug('The given subset was not valid due to quadtratic coupling '
                         'and/or QSS species as specific collision partner')
            logger.debug('This coupling-free subset:' + str(
                species_qss_names_proposition) + 'has been set to replace the previous one')

    return new_ctmech, species_qss_names_proposition
