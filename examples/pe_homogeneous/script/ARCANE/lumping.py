"""Module implementing lumping of species"""

import ARCANE.ct2cti as ct2cti
import ARCANE.tools as tools
import ARCANE.display as display
import cantera as ct
import numpy as np
import scipy
import networkx as nx

logger = display.Logger()
logger.set_log('logCompute')


def find_isomers(mechanism):
    """Function finding isomer species

    Parameters
    ----------
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    isomers_dict_clone : dict
        dictionary containing the isomer species

    """

    isomers_dict = dict()

    species = mechanism.ctmech.species()

    for index_spec, spec in enumerate(species):

        composition = repr(spec.composition)

        if composition not in isomers_dict:
            isomers_dict[composition] = [spec.name]

        else:
            isomers_dict[composition].append(spec.name)

    isomers_dict_clone = isomers_dict.copy()

    # Cleaning species with no isomer
    for key in isomers_dict:

        value = isomers_dict[key]

        if len(value) == 1:
            isomers_dict_clone.pop(key)

    return isomers_dict_clone


def find_weights(cases_list, mechanism, species_names):
    """Calculating the weights of a set of species as ratio of the maxima

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class `ARCANE.mechanisms.Mechanism` object
    species_names :
        list of species

    Returns
    -------
    species_weights : list
        list of species weights

    """

    # Sorting species according to the ctmech
    species_names = tools.sort_species_by_solution(species_names, mechanism.ctmech)

    max_values = []

    for case in cases_list:

        data_dict = case.data_dict(mechanism)

        case_max_values = []
        sum_max = 0

        for spec in species_names:
            profile = data_dict[spec]
            sum_value = sum(profile)
            case_max_values.append(sum_value)
            sum_max += sum_value

        # Normalizing the extracted maximum values
        case_max_values = [value / sum_max for value in case_max_values]

        # Appending to the global list
        max_values.append(case_max_values)

    # Taking the mean value of the weights
    mean_values = np.mean(max_values, 0)

    # Final normalizaition
    sum_means = sum(mean_values)
    species_weights = [value / sum_means for value in mean_values]

    return species_weights


def find_species_ratios(cases_list, mechanism, species_names):
    """Calculating the ratio between species in the cases

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class `ARCANE.mechanisms.Mechanism` object
    species_names :
        list of species

    Returns
    -------
    temp_list : list
        temperature list corresponding to the cases
    spec_ratios_dict : dict
        dictionary containing species ratios

    """

    # Sorting species according to the ctmech
    species_names = tools.sort_species_by_solution(species_names, mechanism.ctmech)

    temp_list = []

    spec_ratios_dict = {}
    for spec in species_names:
        spec_ratios_dict[spec] = []

    for case in cases_list:

        # Extracting the data
        data_dict = case.data_dict(mechanism, reload=True)

        # Loop on every grid point
        for index, grid_point in enumerate(data_dict['Grid']):

            temp = data_dict['T'][index]

            # Store only if not present already
            if temp not in temp_list:

                sum_values = 0
                store_values = []

                for spec in species_names:
                    spec_mass_fraction = data_dict[spec][index]
                    if spec_mass_fraction < 1.0e-15:
                        spec_mass_fraction = 0

                    store_values.append(spec_mass_fraction)

                    sum_values += spec_mass_fraction

                if sum_values != 0:
                    # Normalizing the extracted maximum values
                    store_values = [value / sum_values if value != 0 else 1e-15 for value in store_values]
                    temp_list.append(temp)

                    for index_spec, spec in enumerate(species_names):
                        spec_ratios_dict[spec].append(store_values[index_spec])

    return temp_list, spec_ratios_dict


def automatic_lumped_species_name(mechanism, species_to_lump_names):
    """Creates a name for the species based on its composition and the previous species with the same name
    C*H*O*N*(L*)

    Parameters
    ----------
    mechanism :
        class `ARCANE.mechanisms.Mechanism` object
    species_to_lump_names :
        names of species to be lumped

    Returns
    -------
    lumped_species_name_addition : str
        new species name

    """

    ctmech = mechanism.ctmech
    species = ctmech.species()
    species_names = ctmech.species_names

    species_to_lump = [spec for spec in species if spec.name in species_to_lump_names]

    composition = species_to_lump[0].composition

    lumped_species_name = ''

    if 'C' in composition:
        if composition['C'] > 1:
            lumped_species_name += 'C' + str(int(composition['C']))
        else:
            lumped_species_name += 'C'
    if 'H' in composition:
        if composition['H'] > 1:
            lumped_species_name += 'H' + str(int(composition['H']))
        else:
            lumped_species_name += 'H'
    if 'O' in composition:
        if composition['O'] > 1:
            lumped_species_name += 'O' + str(int(composition['O']))
        else:
            lumped_species_name += 'O'
    if 'N' in composition:
        if composition['N'] > 1:
            lumped_species_name += 'N' + str(int(composition['N']))
        else:
            lumped_species_name += 'N'

    # Verifying if an other lumped species with the same name exists or not
    index = 1

    # Finding correct index
    existent_species = [spec_name for spec_name in species_names
                        if spec_name.startswith(lumped_species_name + '(L')]

    if existent_species:
        lumping_indices = [int(spec.replace(lumped_species_name + '(L', '').replace(')', '')) for spec in
                           existent_species]
        index = max(lumping_indices) + 1

    lumped_species_name_addition = lumped_species_name + '(L' + str(index) + ')'

    return lumped_species_name_addition


def lump_species(mechanism, species_arrhenius_dict, lumped_species_name='',
                 lumped_cti_file='', optimize_NASA=True):
    """Lumps the set of species making a new mechanism object

    Parameters
    ----------
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object
    species_arrhenius_dict :
        dictionary of species Arrhenius ratios
    lumped_species_name :
        specific name for lumped species (if nothing it will be generated automatically) (Default value = '')
    lumped_cti_file :
        name of the output cti file (not written as default)
    optimize_NASA :
        if False, the NASA of the first species will be taken (Default value = True)

    Returns
    -------
    new_ctmech : class `Cantera.Solution` object
        Cantera Solution object with species lumping

    """

    ctmech = mechanism.ctmech
    # Checking from transport data
    if mechanism.transport == 'Transport':
        build_transport = False
    else:
        build_transport = True
    species = ctmech.species()

    # Sorting species according to the ctmech
    species_to_lump_names = tools.sort_species_by_solution(list(species_arrhenius_dict.keys()), mechanism.ctmech)

    species_to_lump = [spec for spec in species if spec.name in species_to_lump_names]

    # If single value spec_arrhenius_dict, transform in Arrhenius format
    if type(species_arrhenius_dict[species_to_lump_names[0]]) in [float, int]:
        sum_values = sum(list(species_arrhenius_dict.values()))
        for spec in species_arrhenius_dict:
            species_arrhenius_dict[spec] = [species_arrhenius_dict[spec] / sum_values, 0, 0]

    # Name of lumped species
    if not lumped_species_name:
        lumped_species_name = automatic_lumped_species_name(mechanism, species_to_lump_names)

    # Create a copy of the first species to lump
    if optimize_NASA:
        high_7_coefficients, low_7_coefficients = tools.fit_NASA(mechanism, spec_arrhenius_dict=species_arrhenius_dict)
        NASA_polynomials = create_NASA_object(300, 1000, 5000, 1.01325e5, low_7_coefficients, high_7_coefficients)
    else:
        NASA_polynomials = None

    lumped_species = create_species(lumped_species_name, species_to_lump[0], NASA=NASA_polynomials,
                                    transport=build_transport)

    new_species = species
    new_species.append(lumped_species)

    # Modifying the reactions
    new_ctmech = ctmech
    for index_spec, spec_name in enumerate(species_to_lump_names):
        arrhenius_parameters = species_arrhenius_dict[spec_name]

        new_reactions = tools.change_reaction_by_species(new_ctmech, spec_name,
                                                         A=arrhenius_parameters[0],
                                                         b=arrhenius_parameters[1],
                                                         Ea=arrhenius_parameters[2],
                                                         new_names_list=lumped_species_name)

        new_species.remove(species_to_lump[index_spec])

        new_ctmech = ct.Solution(thermo='IdealGas',
                                 kinetics='GasKinetics',
                                 species=new_species,
                                 reactions=new_reactions)

    new_ctmech = sum_similar_reactions(new_ctmech)

    # Writes the cti file if specified
    if lumped_cti_file:
        ct2cti.write(new_ctmech, lumped_cti_file)

    return new_ctmech


def create_NASA_object(T_low, T_mid, T_high, P_ref, low_7_coefficients, high_7_coefficients):
    """Creates a NASA polynomial object

    Parameters
    ----------
    T_low :
        low temperature threshold for low coefficients
    T_mid :
        mid temperature threshold common to the to coefficients sets
    T_high :
        high temperature threshold for high coefficients
    P_ref :
        reference pressure
    low_7_coefficients :
        list of 7 coefficients for the low temperature part
    high_7_coefficients :
        list of 7 coefficients for the high temperature part

    Returns
    -------
    NASA : class `Cantera.NasaPoly2` object
        NasaPoly2 Cantera object

    """

    coeffs = list()
    coeffs.append(T_mid)
    coeffs += high_7_coefficients
    coeffs += low_7_coefficients

    NASA = ct.NasaPoly2(T_low, T_high, P_ref, coeffs)

    return NASA


def create_species(name, initial_species, NASA=None, transport=True):
    """Creating a new species from an existing one

    Parameters
    ----------
    name :
        name of the new species
    initial_species :
        base species to copy properties
    NASA :
        user specified NASA polynomial ctmech.species()[0] object (Default value = None)
    transport :
        if False, does not build the transport data (Default value = True)

    Returns
    -------
    new_species : class `Cantera.Species` object
        new class `Cantera.Species` object with the given thermodynamic and transport data

    """

    new_species = ct.Species(name, initial_species.composition,
                             initial_species.charge,
                             initial_species.size)

    if transport:
        new_species.transport = initial_species.transport

    if NASA is None:
        new_species.thermo = initial_species.thermo
    else:
        new_species.thermo = NASA

    return new_species

def lump_non_isomer_species(mechanism, ratios_dict, lumped_species_name=''):
    """Lumps the set of species making a new species

    Parameters
    ----------
    mechanism :
        class `ARCANE.mechanisms.Mechanism` object
    ratios_dict :
        dictionary of species ratios
    lumped_species_name :
        specific name for lumped species (if nothing it will be generated automatically) (Default value = '')

    Returns
    -------
    new_species : class `Cantera.Species` object
        new class `Cantera.Species` object with the given thermodynamic and transport data

    """

    ctmech = mechanism.ctmech
    # Checking from transport data
    if mechanism.transport == 'Transport':
        build_transport = False
    else:
        build_transport = True

    species = ctmech.species()

    # Sorting species according to the ctmech
    species_to_lump_names = tools.sort_species_by_solution(list(ratios_dict.keys()), mechanism.ctmech)

    species_to_lump = [spec for spec in species if spec.name in species_to_lump_names]

    # Name of lumped species
    if not lumped_species_name:
        lumped_species_name = automatic_lumped_species_name(mechanism, species_to_lump_names)

    # Create a copy of the first species to lump
    high_7_coefficients, low_7_coefficients = tools.fit_NASA(mechanism, composition_dict=ratios_dict)
    NASA_polynomials = create_NASA_object(300, 1000, 5000, 1.01325e5, low_7_coefficients, high_7_coefficients)

    new_composition = {'C': 0, 'H': 0, 'O': 0}
    for spec in ratios_dict:
        new_composition['C'] += ctmech.n_atoms(spec, 'C') * ratios_dict[spec]
        new_composition['H'] += ctmech.n_atoms(spec, 'H') * ratios_dict[spec]
        new_composition['O'] += ctmech.n_atoms(spec, 'O') * ratios_dict[spec]

    charge = max([spec.charge for spec in species_to_lump])
    size = max([spec.size for spec in species_to_lump])

    new_species = ct.Species(lumped_species_name, new_composition, charge, size)

    if build_transport:
        new_species.transport = species_to_lump[0].transport

    new_species.thermo = NASA_polynomials

    return new_species

def sum_similar_reactions(ctmech):
    """Adding reactions that have the same equation, same temperature exponent and same activation energy

    Parameters
    ----------
    ctmech :
        class `Cantera.Solution` object


    Returns
    -------
    new_ctmech : class `Cantera.Solution` object
        new class `Cantera.Solution` object with the modifications


    """

    reactions = [tools.copy_reaction(reac) for reac in ctmech.reactions()]

    new_reactions = []
    summed_reactions = []
    duplicate_reactions = []

    for i, reaction_i in enumerate(reactions):
        if i not in summed_reactions:

            for j, reaction_j in enumerate(reactions[i+1:]):
                j += 1

                if reaction_i.equation == reaction_j.equation:

                    if reaction_i.reaction_type in [1, 2] and reaction_j.reaction_type in [1, 2] \
                            and reaction_i.rate.temperature_exponent == reaction_j.rate.temperature_exponent \
                            and reaction_i.rate.activation_energy == reaction_j.rate.activation_energy:
                        reaction_i = tools.change_rate(reaction_i, A_addition=reaction_j.rate.pre_exponential_factor)
                        reaction_i.duplicate = False

                        summed_reactions.append(j + i)

                    else:
                        duplicate_reactions += [i, j + i]

                elif (reaction_i.reactants == reaction_j.reactants
                      and reaction_j.products == reaction_i.products) \
                        or (reaction_i.reactants == reaction_j.products
                            and reaction_j.reactants == reaction_i.products):

                    duplicate_reactions += [i, j + i]

        if i not in summed_reactions:
            if i in duplicate_reactions:
                reaction_i.duplicate = True

            new_reactions.append(reaction_i)

        reactions[i] = reaction_i

    # Creates new Solution object
    new_ctmech = ct.Solution(thermo='IdealGas',
                             kinetics='GasKinetics',
                             species=ctmech.species(),
                             reactions=new_reactions)

    return new_ctmech


def detect_identical_species(mechanism, species_to_lump_names):
    """Detects species with same NASA polynomials inside a list

    Parameters
    ----------
    mechanism :
        class `ARCANE.mechanisms.Mechanism` object
    species_to_lump_names :
        list of species that will be lumped

    Returns
    -------
    groups_to_lump : list
        list of grouped identical species

    """

    ctmech = mechanism.ctmech

    used_species = []

    groups_to_lump = []

    # Loop to gather the NASA coefficients of the first species

    for i, spec_i in enumerate(species_to_lump_names):

        species_to_lump_list = []
        # species_to_lump_list.append(species_to_lump_names[i])

        if spec_i not in used_species:
            used_species.append(species_to_lump_names[i])
            NASA_coeffs_species = ctmech.species(spec_i).thermo.coeffs

            # Loop to gather the NASA coefficients of the second species

            for j, spec_j in enumerate(species_to_lump_names[i + 1:]):

                if spec_j not in used_species:

                    NASA_coeffs_species_tested = ctmech.species(spec_j).thermo.coeffs

                    # Comparison of NASA coefficients

                    if all(NASA_coeffs_species == NASA_coeffs_species_tested):
                        if spec_i not in species_to_lump_list:
                            species_to_lump_list.append(spec_i)
                        if spec_j not in species_to_lump_list:
                            species_to_lump_list.append(spec_j)
                        used_species.append(spec_j)

            if species_to_lump_list:
                groups_to_lump.append(species_to_lump_list)

    return groups_to_lump


def lump_species_2_by_2(mechanism, species_arrhenius_dict, max_error):
    """Function determining the lumping potential of species through thermodynamic proximity

    Parameters
    ----------
    mechanism :
        class `ARCANE.Mecanism` object
    species_arrhenius_dict :
        dictionary of species Arrhenius ratios
    max_error :
        percentage of error authorized between the curves of Cp and Cp_fit
        and enthalpy and enthalpy_fit for lumping

    Returns
    -------
    lumpable_species_2_by_2 : list
        list of list of species to lump together

    """

    ctmech = mechanism.ctmech
    species = ctmech.species()

    # Sorting species according to the ctmech
    species_to_lump_names = tools.sort_species_by_solution(list(species_arrhenius_dict.keys()), mechanism.ctmech)

    # If single value spec_arrhenius_dict, transform in Arrhenius format
    if type(species_arrhenius_dict[species_to_lump_names[0]]) in [float, int]:
        sum_values = sum(list(species_arrhenius_dict.values()))
        for spec in species_arrhenius_dict:
            species_arrhenius_dict[spec] = [species_arrhenius_dict[spec] / sum_values, 0, 0]

    Tmin = 300.
    Tmax = 5000.
    T_rupture = 1000
    temperature_low = np.linspace(Tmin, T_rupture, 201)
    temperature_high = np.linspace(T_rupture, Tmax, 201)

    lumpable_species_2_by_2 = []
    used_species = []

    for i, spec_i in enumerate(species_to_lump_names):

        if species_to_lump_names[i] not in used_species:

            used_species.append(species_to_lump_names[i])

            # Gathering of Cp, enthalpy and entropy of species i alone

            NASA_coeffs_species = ctmech.species(spec_i).thermo.coeffs

            cp_1_low = tools.NASA(temperature_low, *NASA_coeffs_species[1:8], 'cp')
            enthalpy_1_low = tools.NASA(temperature_low, *NASA_coeffs_species[1:8], 'enthalpy')
            entropy_1_low = tools.NASA(temperature_low, *NASA_coeffs_species[1:8], 'entropy')

            cp_1_high = tools.NASA(temperature_high, *NASA_coeffs_species[8:], 'cp')
            enthalpy_1_high = tools.NASA(temperature_high, *NASA_coeffs_species[8:], 'enthalpy')
            entropy_1_high = tools.NASA(temperature_high, *NASA_coeffs_species[8:], 'entropy')

            for j, spec_j in enumerate(species_to_lump_names[i + 1:]):
                j = j + i + 1

                if i != j:

                    species_to_lump_list = []

                    NASA_coeffs_species_tested = ctmech.species(spec_j).thermo.coeffs

                    red_species_arrhenius_dict = {}
                    red_species_arrhenius_dict[spec_i] = species_arrhenius_dict[spec_i]
                    red_species_arrhenius_dict[spec_j] = species_arrhenius_dict[spec_j]

                    # Gathering of Cp, enthalpy and entropy of species j alone

                    cp_2_low = tools.NASA(temperature_low, *NASA_coeffs_species_tested[1:8], 'cp')
                    enthalpy_2_low = tools.NASA(temperature_low, *NASA_coeffs_species_tested[1:8], 'enthalpy')
                    entropy_2_low = tools.NASA(temperature_low, *NASA_coeffs_species_tested[1:8], 'entropy')

                    cp_2_high = tools.NASA(temperature_high, *NASA_coeffs_species_tested[8:], 'cp')
                    enthalpy_2_high = tools.NASA(temperature_high, *NASA_coeffs_species_tested[8:], 'enthalpy')
                    entropy_2_high = tools.NASA(temperature_high, *NASA_coeffs_species_tested[8:], 'entropy')

                    # Calculation of Cp, enthalpy and entropy of species for a mixture of species i and j

                    coeff_low, coeff_high = tools.fit_NASA(mechanism, spec_arrhenius_dict=red_species_arrhenius_dict)

                    cp_fit_low = tools.NASA(temperature_low, *coeff_low, 'cp')
                    enthalpy_fit_low = tools.NASA(temperature_low, *coeff_low, 'enthalpy')
                    entropy_fit_low = tools.NASA(temperature_low, *coeff_low, 'entropy')

                    cp_fit_high = tools.NASA(temperature_high, *coeff_high, 'cp')
                    enthalpy_fit_high = tools.NASA(temperature_high, *coeff_high, 'enthalpy')
                    entropy_fit_high = tools.NASA(temperature_high, *coeff_high, 'entropy')

                    # If error between the mixture Cp to individual Cp,
                    # idem for enthalpy < 10% then species i and j can be lumped together
                    values_to_test = [zip(cp_1_low, cp_fit_low),
                                      zip(enthalpy_1_low, enthalpy_fit_low),
                                      zip(cp_1_high, cp_fit_high),
                                      zip(enthalpy_1_high, enthalpy_fit_high),
                                      zip(cp_2_low, cp_fit_low),
                                      zip(enthalpy_2_low, enthalpy_fit_low),
                                      zip(cp_2_high, cp_fit_high),
                                      zip(enthalpy_2_high, enthalpy_fit_high)]

                    valid = True
                    for set_to_test in values_to_test:
                        error = max([abs((ref - new) / max(abs(ref), abs(new))) for ref, new in set_to_test])

                        if error > max_error:
                            valid = False

                    if valid:
                        species_to_lump_list.append(species_to_lump_names[i])
                        species_to_lump_list.append(species_to_lump_names[j])

                        lumpable_species_2_by_2.append(species_to_lump_list)

    return lumpable_species_2_by_2


def get_all_recombinations(graph):
    """Graph analysis of combination clustering

    Parameters
    ----------
    graph :
        Graph of the species to get the combinations possible

    """
    l = [set(clique) for clique in nx.enumerate_all_cliques(graph)]

    for i in reversed(range(len(l))):

        current_species_used = l[i]
        temp_res = [sorted(list(l[i]))]

        for j in reversed(range(len(l))):

            if l[j].isdisjoint(current_species_used):
                current_species_used |= l[j]
                temp_res.append(list(sorted(l[j])))

        if len(current_species_used) == len(graph.nodes):
            yield temp_res


def get_best_combination(lumpable_species_2_by_2):
    """Selects the best clustering in a list of species couples

    Parameters
    ----------
    lumpable_species_2_by_2 :
        list of lumpable species 2 by 2

    Returns
    -------
    best_combination : list
        list of list of possible species best combinations for lumping

    """

    G = nx.Graph()
    best_combination = []

    G.add_edges_from(lumpable_species_2_by_2)

    result = sorted(get_all_recombinations(G), key=len)
    min_len = len(result[0])

    i = 0

    while len(result[i]) == min_len:
        i += 1
        best_combination.append(result[i - 1])

    return best_combination


def fit_sum_arrhenius_reac(temperature_list, list_of_coefficients):
    """Fits the rates constants sum of several reaction with a single set of Arrhenius coefficients

    Parameters
    ----------
    temperature_list :
        list of temperature to perform the fit on
    list_of_coefficients :
        list of all the Arrhenius coefficients to fit

    Returns
    -------
    new_coefficients : list
        new coefficients of the Arrhenius

    """

    k = 0

    for i in range(len(list_of_coefficients)):
        k += tools.Arrhenius_function(temperature_list, *list_of_coefficients[i])

    new_coefficients = tools.Arrhenius_fit(temperature_list, k)

    return new_coefficients


def create_spec_arrhenius_dict(temp_list, spec_ratio_dict, option='best'):
    """Creates a dictionary with species as keys and Arrhenius coefficient representing the composition evolution as values

    Parameters
    ----------
    temp_list :
        list of temperatures
    spec_ratio_dict :
        dictionary with species as keys and ratios as values
    option :
        specifies which type of function will be used (Arrhenius, Arrhenius with T dependence or constant) (Default value = 'best')

    Returns
    -------
    spec_arrhenius_dict : dict
        dictionary of species Arrhenius coefficients

    """
    spec_arrhenius_dict = {}
    for spec in spec_ratio_dict:
        logger.debug('Fitting Arrhenius for species ' + spec)
        arrh = tools.Arrhenius_fit(temp_list, spec_ratio_dict[spec], option=option)
        if arrh[0] < 0:
            arrh[0] = 1
        spec_arrhenius_dict[spec] = arrh

    return spec_arrhenius_dict
