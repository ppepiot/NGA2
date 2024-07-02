"""Miscellaneous tools and utilities"""

# Import statements
import os
import re
import shutil
import sys

import cantera as ct
import numpy as np
import scipy.optimize as optimize
from scipy.integrate import simps
from scipy.signal import find_peaks

import ARCANE.cases as cases
import ARCANE.display as display
import ARCANE.kwdict as kwdict
import ARCANE.sampling as sampling
import ARCANE.database as database

kwdict = kwdict.Kwdict()

logger = display.Logger()
logger.set_log('logCompute')


def set_composition_from_string(text, ctmech, dict_in_mass=False):
    """Convert text contained in entry into a mole-fraction dictionary.
    Types allowed are X and Y. Numerical values are normalized automatically.
    All species must belong to mechanism.

    Parameters
    ----------
    text : str
        composition, formatted as ``'type/species_1/amount_1/species_2/amount_2...'``
    ctmech :`Cantera.Solution`
        class `Cantera.Solution` object
    dict_in_mass : bool
        (Default value = False)

    Returns
    -------
     : dict
        dictionary containing the mole fraction for different species

    Example
    -------
    >>> set_composition_from_string('X/CH4/1', ctmech)
    >>> set_composition_from_string('X/O2/0.21/N2/0.79', ctmech)

    """

    if type(text) == dict:
        composition = text

        if dict_in_mass:
            ctmech.Y = composition
            composition = ctmech.mole_fraction_dict()
            
    elif text in ctmech.species_names:
        composition = {text: 1}
    else:
        try:
            # Convert text to list
            entry = text.split('/')

            # Type of data is first parameter: if first entry is not X and Y, raise an exception
            mytype = entry.pop(0)

            # Species list (even entries): if not in species list, raise an exception
            slist = entry[0::2]
            # Amounts list (odd entries): if not positive numerical values, raise an exception
            xlist = entry[1::2]
            # Convert xlist to floats and normalize numerical values
            xlist = [float(myf) for myf in xlist]
            xlist = [myf / sum(xlist) for myf in xlist]

        except AssertionError:
            logger.error(
                "Composition entries must be formatted as type_of_amount/Species_1/Amount_1/Species_2/Amount_2...")

        try:
            # Check that all species are in mechanism
            s1 = set(ctmech.species_names)
            s2 = set(slist)
            assert s2.issubset(s1)

        except AssertionError:
            logger.error("Species in composition {0} are not part of the mechanism".format(text))

        # If type is 'Y', convert to mole fractions using cantera
        if mytype == 'Y':
            ctmech.Y = dict(zip(slist, xlist))
            composition = ctmech.mole_fraction_dict()
        else:
            composition = dict(zip(slist, xlist))

    # Return dictionary, skipping 0
    sum_of_of_values = sum(composition.values())

    return {k: round(v, 15) / sum_of_of_values for k, v in composition.items() if v}


def set_composition_from_phi(phi, fuel, oxidizer, ctmech):
    """Convert phi into mole-fraction dictionary according to specified fuel and oxidizer dicts.

    Parameters
    ----------
    phi : float
        equivalence ratio :math:`\phi`
    fuel : dict
        dict containing fuel definition
    oxidizer : dict
        dict containing oxidizer definition
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object

    Returns
    -------
     : dict
        dictionary containing the mole fraction for different species

    TODO
    ----
    Write unittest for function ``set_composition_from_phi()`` and test function

    """

    soxi = compute_stoichiometric_ratio(fuel, oxidizer, ctmech, type="mole")

    # Initialize
    X = {}

    effective_oxidizer = [key for key in oxidizer.keys() if key not in ['N2']]
    effective_oxidizer = effective_oxidizer[0]

    X[effective_oxidizer] = oxidizer[effective_oxidizer] * soxi

    # Mole fractions of oxidizer
    for key in oxidizer:
        X[key] = oxidizer[key] * X[effective_oxidizer]

    # Mole fractions of fuel
    for key in fuel:
        if key in X:
            X[key] += fuel[key] * phi / soxi * X[effective_oxidizer]
        else:
            X[key] = fuel[key] * phi / soxi * X[effective_oxidizer]

    composition = X

    # Return dictionary, skipping 0
    return {k: round(v, 8) if v > 0 else 0 for k, v in composition.items()}


def set_phi_from_composition(X, fuel, oxidizer, ctmech):
    """Compute equivalence ratio :math:`\phi` corresponding to mole-fraction dictionary
    according to specified fuel and oxidizer dicts.

    Parameters
    ----------
    X : dict
        dict containing composition of interest in mole fraction
    fuel : dict
        dict containing fuel definition
    oxidizer : dict
        dict containing oxidizer definition
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object


    Returns
    -------
    soxi / X_oxi : float
        equivalence ratio :math:`\phi` corresponding to the given composition

    TODO
    ----
    - Write unittest for function ``set_phi_from_composition()`` and test function
    - Handle exceptions more properly


    """

    # Get total moles of fuels
    XF = 0.0
    for key in fuel:
        if not X[key]:
            text = "Fuel species " + str(key) + " not in provided composition " + str(X)
            sys.exit(text)
        XF += X[key]

    # Make sure composition is compatible with fuel definition
    for key in fuel:
        if abs(X[key] / XF - fuel[key] / sum(fuel.values())) > 1e-2:
            text = "Fuel definition " + key + ' : ' + str(fuel[key]) \
                   + " does not match provided composition " + str(X)
            sys.exit(text)

    # Get total moles of oxidizer
    XO = 0.0
    for key in oxidizer:
        if not X[key]:
            text = "Oxidizer species " + str(oxidizer) + " not in provided composition " + str(X)
            sys.exit(text)
        XO += X[key]

    # Make sure composition is compatible with fuel definition
    for key in oxidizer:
        if abs(X[key] / XO - oxidizer[key] / sum(oxidizer.values())) > 1e-2:
            text = "Oxidizer definition " + key + ' : ' + str(oxidizer[key]) \
                   + " does not match provided composition " + str(X)
            sys.exit(text)

    effective_oxidizer = [key for key in oxidizer.keys() if key not in ['N2']]
    effective_oxidizer = effective_oxidizer[0]

    oxi_ratio = 1 / oxidizer[effective_oxidizer]

    for key in oxidizer:
        oxidizer[key] = oxidizer[key] * oxi_ratio

    # Creating a dictionary with composition
    composition_dict = {'C': {}, 'H': {}, 'Al': {}, 'O': {}}
    elements = ctmech.element_names
    for element in composition_dict:
        for species in ctmech.species_names:
            if element in elements:
                composition_dict[element][species] = ctmech.n_atoms(species, element)
            else:
                composition_dict[element][species] = 0

    # Compute phi for oxidizer of type CxHyOz
    soxi_num = 0.0
    for key in fuel:
        soxi_num += X[key] * (- 4 * composition_dict['C'][key]
                              - composition_dict['H'][key]
                              - 3 * composition_dict['Al'][key]
                              + 2 * composition_dict['O'][key])

    soxi_denom = 0.0
    for key in oxidizer:
        soxi_denom += oxidizer[key] * (4 * composition_dict['C'][key]
                                       + composition_dict['H'][key]
                                       + 3 * composition_dict['Al'][key]
                                       - 2 * composition_dict['O'][key])

    soxi = soxi_num / soxi_denom

    if soxi < 0:
        soxi = 1

    X_oxi = X[effective_oxidizer]

    return soxi / X_oxi


def set_composition_from_new_phi(phi, fuel, oxidizer, ctmech):
    """
    Convert phi into mole-fraction dictionary according to specified fuel and oxidizer dicts.
    This function uses an atomic conservation formula of the equivalence ratio equivalent to the classical one
    for CxHy hydrocarbons but yielding different results for oxygenated fuels CxHyOz.

    Parameters
    ----------
    phi : float
        equivalence ratio :math:`\phi`
    fuel : dict
        dict containing fuel definition
    oxidizer : dict
        dict containing oxidizer definition
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object

    Returns
    -------
     : dict
        dictionary containing the mole fraction for different species

    TODO
    ----
    Write unittest for function ``set_composition_from_phi()`` and test function

    """
    elements_in_mech = ctmech.element_names
    composition_dict = {element: {} for element in elements_in_mech}
    for element in composition_dict:
        for species in ctmech.species_names:
            composition_dict[element][species] = ctmech.n_atoms(species, element)

    # Normalizing fuel
    sum_of_fuels = sum(fuel.values())
    for key in fuel:
        fuel[key] = fuel[key] / sum_of_fuels

    # Normalizing oxidizer
    sum_of_fuels = sum(oxidizer.values())
    for key in oxidizer:
        oxidizer[key] = oxidizer[key] / sum_of_fuels

    # Compute phi
    fuel_atom_dict = {kwdict.elements_names[element]: 0 for element in elements_in_mech}
    oxidizer_atom_dict = {kwdict.elements_names[element]: 0 for element in elements_in_mech}

    phi_numerator = 0
    phi_denominator = 0
    oxidizer_coefficient = 0
    oxidizer_in_fuel = False
    fuel_in_oxidizer = False

    # Determination of the fuel composition
    for key in fuel:
        for element in elements_in_mech:
            if element not in kwdict.elements_names:
                kwdict.elements_names[element] = 'unknown element'
                kwdict.elements_weights = 0
                kwdict.elements_phi_coeff = 0

            fuel_atom_dict[kwdict.elements_names[element]] += fuel[key] * composition_dict[element][key]

            if not oxidizer_in_fuel and kwdict.elements_phi_coeff[element] < 0 \
                    and composition_dict[element][key] > 0:
                oxidizer_in_fuel = True

    # Determination of the oxidizer composition
    for key in oxidizer:
        for element in elements_in_mech:
            if element not in kwdict.elements_names:
                kwdict.elements_names[element] = 'unknown element'
                kwdict.elements_weights = 0
                kwdict.elements_phi_coeff = 0

            oxidizer_atom_dict[kwdict.elements_names[element]] += composition_dict[element][key]

            if not fuel_in_oxidizer and kwdict.elements_phi_coeff[element] > 0 \
                    and composition_dict[element][key] > 0:
                fuel_in_oxidizer = True
                reducers = [element for element in kwdict.elements_phi_coeff if kwdict.elements_phi_coeff[element] > 0]
                logger.warning(f"A reducer ({','.join(reducers)} is part of your oxidizer "
                               f"and this is not yet supported.")

    for element in elements_in_mech:
        if kwdict.elements_phi_coeff[element] > 0:
            phi_numerator += kwdict.elements_phi_coeff[element] * fuel_atom_dict[f'{kwdict.elements_names[element]}']
        else:
            phi_denominator -= kwdict.elements_phi_coeff[element] * fuel_atom_dict[f'{kwdict.elements_names[element]}']
            oxidizer_coefficient -= kwdict.elements_phi_coeff[element] \
                                    * oxidizer_atom_dict[f'{kwdict.elements_names[element]}']

    if oxidizer_in_fuel:
        max_phi = phi_numerator / phi_denominator
        if phi > max_phi:
            phi = max_phi
            logger.debug(f'For an oxygenated fuel, the equivalence ratio is bounded by a maximum value defined as '
                           f'the equivalence ratio of the sole fuel, {max_phi} is the present case.\n'
                           f'The composition will be set with this maximum value (i.e. without oxidizer).')

    oxidizer_ratio = (1 / oxidizer_coefficient) * ((phi_numerator / phi) - phi_denominator)

    # Initialize
    X = {}

    effective_oxidizer = [key for key in oxidizer.keys() if key not in ['N2']]
    effective_oxidizer = effective_oxidizer[0]

    # Mole fractions of oxidizer
    for key in oxidizer:
        X[key] = oxidizer[key] / oxidizer[effective_oxidizer] * oxidizer_ratio

    # Mole fractions of fuel
    sum_of_fuels = sum(fuel.values())
    for key in fuel:
        X[key] = fuel[key] / sum_of_fuels

    composition = X

    # Return dictionary, skipping 0
    return {k: round(v, 8) if v > 0 else 0 for k, v in composition.items()}


def set_new_phi_from_composition(X, ctmech):
    """
    Compute equivalence ratio :math:`\phi` corresponding to mole-fraction dictionary
    according to specified fuel and oxidizer dicts.
    This function uses an atomic conservation formula of the equivalence ratio equivalent to the classical one
    for CxHy hydrocarbons but yielding different results for oxygenated fuels CxHyOz.

    Parameters
    ----------
    X : dict
        dict containing composition of interest in mole fraction
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object


    Returns
    -------
    phi : float
        equivalence ratio :math:`\phi` corresponding to the given composition

    TODO
    ----
    - Write unittest for function ``set_phi_from_composition()`` and test function
    - Handle exceptions more properly


    """

    elements_in_mech = ctmech.element_names
    composition_dict = {element: {} for element in elements_in_mech}
    for element in composition_dict:
        for species in ctmech.species_names:
            composition_dict[element][species] = ctmech.n_atoms(species, element)

    # Compute phi
    atom_dict = {kwdict.elements_names[element]: 0 for element in elements_in_mech}

    phi_numerator = 0
    phi_denominator = 0

    for key in X:
        for element in elements_in_mech:
            if element not in kwdict.elements_names:
                kwdict.elements_names[element] = 'unknown element'
                kwdict.elements_weights = 0
                kwdict.elements_phi_coeff = 0

            atom_dict[kwdict.elements_names[element]] += X[key] * composition_dict[element][key]

    for element in elements_in_mech:
        if kwdict.elements_phi_coeff[element] > 0:
            phi_numerator += kwdict.elements_phi_coeff[element] * atom_dict[f'{kwdict.elements_names[element]}']
        else:
            phi_denominator -= kwdict.elements_phi_coeff[element] * atom_dict[f'{kwdict.elements_names[element]}']

    phi = phi_numerator / phi_denominator

    return round(phi, 5)


def compute_stoichiometric_ratio(fuel, oxidizer, ctmech, type="mole"):
    """
    
    Parameters
    ----------
    fuel : dict
        dict containing fuel definition
    oxidizer : dict
        dict containing oxidizer definition
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object
    method : str
        Specifies the type of ratio computed: mole, mass or mass_oxygen (Default: "mole")

    Returns
    -------
    stoich_ratio : float
        Specified stoichiometric ratio

    """
    
    effective_oxidizer = [key for key in oxidizer.keys() if key not in ['N2']]
    effective_oxidizer = effective_oxidizer[0]

    oxi_ratio = 1 / oxidizer[effective_oxidizer]

    for key in oxidizer:
        oxidizer[key] = oxidizer[key] * oxi_ratio

    composition_dict = {'C': {}, 'H': {}, 'O': {}, 'Al': {}}
    elements = ctmech.element_names
    for element in composition_dict:
        for species in ctmech.species_names:
            if element in elements:
                composition_dict[element][species] = ctmech.n_atoms(species, element)
            else:
                composition_dict[element][species] = 0

    # Compute phi for oxidizer of type CxHyOz
    stoich_ratio_num = 0.0
    for key in fuel:
        stoich_ratio_num += fuel[key] * (- 4 * composition_dict['C'][key]
                                 - composition_dict['H'][key]
                                 - 3 * composition_dict['Al'][key]
                                 + 2 * composition_dict['O'][key])

    # Compute phi
    stoich_ratio_denom = 0.0
    for key in oxidizer:
        stoich_ratio_denom += oxidizer[key] * (4 * composition_dict['C'][key]
                                       + composition_dict['H'][key]
                                       + 3 * composition_dict['Al'][key]
                                       - 2 * composition_dict['O'][key])

    stoich_ratio = stoich_ratio_num / stoich_ratio_denom

    if type == 'mass':
        oxidizer_molecular_weight = sum([oxidizer[spec] * ctmech.molecular_weights[ctmech.species_names.index(spec)] for spec in oxidizer])
        fuel_molecular_weight = sum([fuel[spec] * ctmech.molecular_weights[ctmech.species_names.index(spec)] for spec in fuel])

        stoich_ratio = stoich_ratio * (oxidizer_molecular_weight / fuel_molecular_weight)

    elif type == 'mass_oxygen':
        oxidizer_molecular_weight = ctmech.molecular_weights[ctmech.species_names.index('O2')]
        fuel_molecular_weight = sum([fuel[spec] * ctmech.molecular_weights[ctmech.species_names.index(spec)] for spec in fuel])

        stoich_ratio = stoich_ratio * (oxidizer_molecular_weight / fuel_molecular_weight)

    if type not in ['mole', 'mass', 'mass_oxygen']:
        logger.error(f"This type of stoichiometric ratio {type} is not valid.")

    return stoich_ratio
    

def set_phi_from_fuel_air_ratio(FAR, fuel, oxidizer, ctmech):
    """Compute equivalence ratio :math:`\phi` corresponding to the fuel_air_ratio
    according to specified fuel and oxidizer dicts.

    Parameters
    ----------
    FAR : float
        value of the fuel-to-air ratio
    fuel : dict
        dict containing fuel definition
    oxidizer : dict
        dict containing oxidizer definition
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object


    Returns
    -------
    phi : float
        equivalence ratio :math:`\phi` corresponding to the given FAR
    
    """
    
    # Computing stoichiometric mass ratio
    mass_stoich_ratio = compute_stoichiometric_ratio(fuel, oxidizer, ctmech, type="mass")

    phi = FAR * mass_stoich_ratio

    return phi


def set_fuel_air_ratio_from_phi(phi, fuel, oxidizer, ctmech):
    """Compute equivalence ratio :math:`\phi` corresponding to the fuel_air_ratio
    according to specified fuel and oxidizer dicts.

    Parameters
    ----------
    phi : float
        equivalence ratio :math:`\phi`
    fuel : dict
        dict containing fuel definition
    oxidizer : dict
        dict containing oxidizer definition
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object


    Returns
    -------
    FAR : float
        value of the fuel-to-air ratio corresponding to the given :math:`\phi`

    """

    # Computing stoichiometric mass ratio
    mass_stoich_ratio = compute_stoichiometric_ratio(fuel, oxidizer, ctmech, type="mass")

    FAR = phi * mass_stoich_ratio

    return FAR


def set_mass_flow_rates(ctmech, phi, strain_rate, fuel, oxidizer, T_fuel, T_oxi, width):
    """Sets fuel and oxidizer mass flow rates according to equivalence ratio :math:`\phi`

    Parameters
    ----------
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object
    phi : flat
        equivalence ratio :math:`\phi`
    strain_rate :
        strain rate :math:`a`
    fuel : dict
        dict containing fuel definition
    oxidizer : dict
        dict containing oxidizer definition
    T_fuel : float
        fuel temperature
    T_oxi : float
        oxidizer temperature
    width : float
        domain width

    Returns
    -------
    mdot_fuel : float
        fuel mass flow rate
    mdot_oxi : float
        oxidizer mass flow rate


    """

    # Discarding diluent species
    effective_oxidizer = [key for key in oxidizer.keys() if key not in ['N2']]
    effective_oxidizer = effective_oxidizer[0]

    oxi_ratio = 1 / oxidizer[effective_oxidizer]

    for key in oxidizer:
        oxidizer[key] = oxidizer[key] * oxi_ratio

    # Compute phi for oxidizer of type CxHyOz
    soxi_num = 0.0
    for key in fuel:
        soxi_num += fuel[key] * (- 4 * ctmech.n_atoms(key, 'C')
                                 - ctmech.n_atoms(key, 'H')
                                 + 2 * ctmech.n_atoms(key, 'O'))

    # Compute phi
    soxi_denom = 0.0
    for key in oxidizer:
        soxi_denom += oxidizer[key] * (4 * ctmech.n_atoms(key, 'C')
                                       + ctmech.n_atoms(key, 'H')
                                       - 2 * ctmech.n_atoms(key, 'O'))

    soxi = soxi_num / soxi_denom

    # Fuel molecular weight
    ctmech.X = fuel
    W_fuel = ctmech.mean_molecular_weight

    # Oxidizer molecular weight
    ctmech.X = oxidizer
    W_oxi = ctmech.mean_molecular_weight

    # Fuel stoechiometric mass flow rate
    mdot_fuel = 1

    # Oxidizer stoechiometric mass flow rate
    mdot_oxi = mdot_fuel * (W_oxi / W_fuel) * soxi

    # Fuel mass flow rate
    mdot_fuel *= phi

    # Fuel inlet velocity
    ctmech.TPX = T_fuel, ctmech.P, fuel
    u_fuel = mdot_fuel / ctmech.density

    # Oxidizer inlet velocity
    ctmech.TPX = T_oxi, ctmech.P, oxidizer
    u_oxi = mdot_oxi / ctmech.density

    # Target strain rate :math:`a`
    target_strain_rate = strain_rate

    # Actual strain rate :math:`a`
    actual_strain_rate = (u_fuel + u_oxi) / width

    # Correction factor for the mass flow rates
    correction = target_strain_rate / actual_strain_rate

    # Correcting the mass flow rates
    mdot_fuel *= correction
    mdot_oxi *= correction

    return mdot_fuel, mdot_oxi


def expand_values(text):
    """Convert short-hand notations to list of values

    Three notations are allowed:

    - ``'-'``: treat as dash-separated list of values
    - ``'/'``: treat as lower bound/higher bound/delta, simple equivalent of 'range' function.

    Parameters
    ----------
    text : str
        input string

    Returns
    -------
    values : list
        list of values corresponding to the text input


    Example
    -------
    >>> expand_values("1000/1500/100")
    [1000, 1100, 1200, 1300, 1400, 1500]
    >>> expand_values("1000-1250-1500")
    [1000, 1250, 1500]
    >>> expand_values("0.8/1.4/0.1-0.95/1.05/0.05")
    [0.8, 0.9, 1.05, 1.0, 0.95, 1.1, 1.2, 1.3, 1.4]


    """

    sort_and_set = False

    # Initialize
    values = []

    if re.search('/', text):
        ranges = text.split('-')
        if len(ranges) > 1:
            sort_and_set = True
        for range in ranges:
            tmp = range.split('/')
            values = np.concatenate([np.arange(
                float(tmp[0]), float(tmp[1]) + float(tmp[2]) / 2, float(tmp[2])), values])
    elif re.search('-', text):
        values = text.split('-')
    else:
        values.append(text)

    # Convert into floats
    try:
        values = [round(float(myv), 10) for myv in values]
    except ValueError:
        logger.error("List of T/P/phi must be formatted according to pre-defined rules. See help for guidance.")

    if sort_and_set:
        values = list(set(np.sort(values)))

    # Return
    return values


def convert_to_valid_fortran_var(mylist):
    """Go through a list of strings and replace fortran-incompatible characters with compatible ones
    Typically used to convert list of species names in a mechanism to valid fortran names

    Parameters
    ----------
    mylist : list
        list of strings to modify

    Returns
    -------
    mylist_valid : list
        list of strings modified in a fortran-friendly style


    """

    # Initialize new list
    mylist_valid = []

    if not isinstance(mylist, list):
        mylist = [mylist]

    for item in mylist:
        # Fortran variables cannot start with a number
        item = re.sub(r'\A([0-9])', r'D\1', item)
        # Replace all - by X
        item = re.sub(r'\-', r'X', item)
        # Fortran variables must only contain 0-9,a-z,A-Z and _, replace all other by G
        item = re.sub(r'(\W)', r'G', item)
        # Switch to all upper case for uniformity
        item = item.upper()
        # Check that it is not over 31 characters
        if len(item) > 31:
            sys.exit('Error in fortran name conversion, variable cannot be over 31 characters:', item)
        mylist_valid.append(item)

    return mylist_valid


def get_ignition_delay_time(case, mechanism=None, method='hr'):
    """Compute ignition delay time from max heat release rate

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : :func:`~ARCANE.mechanisms.Mechanism`
        class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
    method : str
        String to define the method used (Default value = 'hr')


    Returns
    -------
    tig : float
        Ignition delay time


    Created: 17/11/13 [PP]

    Last modified: 18/02/12 [QC]


    """

    if not mechanism:
        mechanism = case.mechanism

    if not method or method == 'none':
        method = 'hr'

    if method in kwdict.names['HR']:

        data_dict = case.data_dict(mechanism)

        axis = data_dict['grid']
        temperature = data_dict['Temperature']
        heat_release = data_dict['HeatRelease']

        init_temperature = temperature[0]
        tig = np.nan

        axis = [value for i, value in enumerate(axis) if heat_release[i] != 0]
        temperature = [value for i, value in enumerate(temperature) if heat_release[i] != 0]
        heat_release = [value for value in heat_release if value != 0]

        if heat_release:
            # Minds all peaks
            maxloc_array = find_peaks(np.array(heat_release))[0]
            if len(maxloc_array) >= 2:
                HR_max = 0
                maxloc = len(heat_release) - 1
                # If several peaks found, keeps the maximum
                for pos in maxloc_array:
                    if heat_release[pos] > HR_max and temperature[pos] > init_temperature + 25:
                        HR_max = heat_release[pos]
                        maxloc = pos
            elif len(maxloc_array) > 0:
                maxloc = maxloc_array[0]
            else:
                maxloc = 0

            if maxloc != len(heat_release) - 1:
                tig = axis[maxloc]

    else:
        logger.error('Your method should be: ' + ', '.join(kwdict.compute_methods['tig']))
        raise NameError

    return tig


def get_laminar_flame_velocity(case, mechanism=None, method='inlet'):
    """Computes the laminar flame velocity of a 1D prexmied flame

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
    method : str
        String to define the method used, might be ``'inlet'``, ``'fuel'``, ``'HR'``, ``'ini'`` or ``'abs'``` (Default value = 'inlet')

    Returns
    -------
    sl : float
        Laminar flame velocity


    """
    if case.reactor_type not in ['C1DP', 'C1DPSOOT'] :
        logger.error('Your case should be a premixed flame to work !')
        return None

    if not mechanism:
        mechanism = case.mechanism

    case.run(mechanism)

    ctmech = mechanism.ctmech

    if not method or method == 'none':
        method = 'inlet'

    if method == 'fuel':
        fuel = list(case.fuel.keys())
        sl = 0

        data_dict = case.data_dict(mechanism)
        rho = data_dict['Density']
        x = data_dict['grid']

        for spec in fuel:
            y_fuel = case.extract_quantity(spec)
            source_term_fuel = np.array(data_dict[spec + ' cdot'])
            molar_mass_fuel = ctmech.molecular_weights[ctmech.species_index(spec)]
            y = source_term_fuel * molar_mass_fuel
            sl += -1 / (rho[0] * (y_fuel[0] - y_fuel[-1])) * simps(source_term_fuel * molar_mass_fuel, x) * case.fuel[
                spec]

    elif method in kwdict.names['HR']:
        hr = case.extract_quantity('HR', mechanism=mechanism)
        x = case.extract_quantity('x', mechanism=mechanism)
        T = case.extract_quantity('T', mechanism=mechanism)
        rho = case.extract_quantity('density', mechanism=mechanism)
        flame_object = case.access_cantera_object()
        cp = flame_object.cp_mass
        sl = 1 / (rho[0] * (T[-1] - T[0])) * simps(hr / cp, x)

    elif method in kwdict.methods['init']:
        sl = case.extract_quantity('u init', mechanism=mechanism)

    elif method in kwdict.methods['abs']:
        rho = case.extract_quantity('density', mechanism=mechanism)
        u = case.extract_quantity('u', mechanism=mechanism)
        sl = u[0] - (rho[1] * u[1] - rho[0] * u[0]) / (rho[1] - rho[0])
    else:
        logger.error('Your method should be: ' + ', '.join(kwdict.compute_methods['sl']))
        raise NameError

    return sl


def get_thickness(case, mechanism=None, method='thermal'):
    """Computes the thickness of a 1D premixed flame

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
    method :
        String to define the method used (Default value = 'thermal')

    Returns
    -------
    thickness : float
        Thermal flame thickness


    """
    if case.reactor_type not in ['C1DP', 'C1DPSOOT'] :
        logger.error('Your case should be a premixed flame to work !')
        quit()

    if not mechanism:
        mechanism = case.mechanism

    if not method or method == 'none':
        method = 'thermal'

    if method == 'thermal':

        x = case.extract_quantity('x', mechanism=mechanism)
        T = case.extract_quantity('T', mechanism=mechanism)
        thickness = (max(T) - min(T)) / max(np.gradient(T, x))

    elif method == 'diffusive':
        thermal_conductivity = case.access_cantera_object().thermal_conductivity
        cp = case.access_cantera_object().cp_mass
        rho = case.extract_quantity('density', mechanism=mechanism)
        thickness = thermal_conductivity[0] / (cp[0] * rho[0]) / get_laminar_flame_velocity(case, mechanism)

    elif method == 'Blint':
        thermal_conductivity = case.access_cantera_object().thermal_conductivity
        cp = case.access_cantera_object().cp_mass
        thickness = 2 * get_thickness(case, mechanism, method="diffusive") \
                    * (thermal_conductivity[-1] / cp[-1]) / (thermal_conductivity[0] / cp[0])

    elif method == 'global':

        x = case.extract_quantity('x', mechanism=mechanism)
        T = case.extract_quantity('T', mechanism=mechanism)

        theta = (T - T[0]) / (T[-1] - T[0])
        x_min_found = False
        x_max_found = False
        x_min = 0
        x_max = 0
        for x_local, theta_local in zip(x, theta):
            if theta_local > 0.01 and not x_min_found:
                x_min = x_local
                x_min_found = True
            if theta_local > 0.99 and not x_max_found:
                x_max = x_local
                x_max_found = True
        if not x_max_found or not x_min_found:
            if not x_min_found:
                logger.error('The minimum has not been found !')
            if not x_max_found:
                logger.error('The maximum has not been found !')
            raise ValueError
        thickness = x_max - x_min

    else:
        logger.error('The keyword should be: ' + ', '.join(kwdict.compute_methods['thickness']))
        raise NameError

    return thickness


def get_species_timescale(cases_list, mechanism, concentration_weighting=False, explicit=False,
                          abs_tol=1e-9, rel_tol=1e-5):
    r"""Computes the residence time of the species.

    In the implicit formulation, the concentration is varied as ``abs_tol + c * rel_tol``

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism : :func:`~ARCANE.mechanisms.Mechanism`
        class :func:`~ARCANE.mechanisms.Mechanism` object
    concentration_weighting : bool
        if True, timescales are multiplied by the species concentration (Default value = False)
    explicit : bool
        if True, uses explicit formulation of the timescales :math:`c / \dot{c}`.
        The explicit method is useful for explicit codes as it gives an idea of the time step necessary
        for a proper tempora resolution.
        The implicit method represent the real time scale of the species
        and is the one used in kinetics analysis, in the selection of QSS species for example. (Default value = False)
    abs_tol : float
        absolute tolerance, if the concentration of the species is inferior to abs_tol.
        its timescale will not be computed (Default value = 1e-9)
    rel_tol : float
        relative tolerance (Default value = 1e-5)

    Returns
    -------
    timescale_dict : dict
        dictionary containing the timescales associated to each species


    """

    # Compiling the mechanism
    if mechanism.f90:
        mechanism.reset()
        mechanism.compile()

    # Creating samples database
    samples_dict = sampling.create_samples_database(cases_list, mechanism, ['T', 'P', 'Y', 'HR'], integrated=False)
    n_data = len(samples_dict['grid'])

    # Shorthand
    ctmech = mechanism.ctmech

    species = list(ctmech.species_names)

    if concentration_weighting:
        min_timescale = np.zeros(len(species))
    else:
        min_timescale = np.ones(len(species))

    # Loop through samples and get residence times
    for index in range(n_data):

        ctmech.TPY = samples_dict['T'][index], samples_dict['P'][index], samples_dict['Y'][index]

        concentration_save = np.zeros(len(ctmech.species()))
        concentration_save[:] = ctmech.concentrations[:]  # concentration

        local_dict = {}
        local_dict['T'] = samples_dict['T'][index]
        local_dict['P'] = samples_dict['P'][index]
        local_dict['Y'] = samples_dict['Y'][index]

        timescale = local_timescale(local_dict, mechanism, concentration_weighting, explicit=explicit,
                                    abs_tol=abs_tol, rel_tol=rel_tol)

        for index_spec, spec in enumerate(species):

            if concentration_weighting:
                # Looking for the maximum of timescale*concentration
                if timescale[index_spec] > min_timescale[index_spec]:
                    min_timescale[index_spec] = timescale[index_spec]
            else:
                # Looking for the minimum timescale
                if timescale[index_spec] < min_timescale[index_spec]:
                    min_timescale[index_spec] = timescale[index_spec]

    timescale_dict = dict(zip(species, min_timescale))

    # Update the mechanism object attribute
    mechanism.timescales = timescale_dict

    return timescale_dict


def local_timescale(samples_dict, mechanism, concentration_weighting=False, explicit=False,
                    abs_tol=1e-9, rel_tol=1e-5):
    r"""Computes the residence time of the species.

    In the implicit formulation, the concentration is varied as ``abs_tol + c * rel_tol``

    Parameters
    ----------
    samples_dict : dict
        data dictionary of sampled data
    mechanism : :func:`~ARCANE.mechanisms.Mechanism`
        class :func:`~ARCANE.mechanisms.Mechanism` object
    concentration_weighting : bool
        if True, timescales are multiplied by the species concentration (Default value = False)
    explicit : bool
        if True, uses explicit formulation of the timescales :math:`c / \dot{c}`.
        The explicit method is useful for explicit codes as it gives an idea of the time step necessary
        for a proper tempora resolution.
        The implicit method represent the real time scale of the species
        and is the one used in kinetics analysis, in the selection of QSS species for example. (Default value = False)
    abs_tol : float
        absolute tolerance, if the concentration of the species is inferior to abs_tol.
        its timescale will not be computed (Default value = 1e-9)
    rel_tol : float
        relative tolerance (Default value = 1e-5)

    Returns
    -------
    timescale : dict
        dictionary of local timescales for each species


    """

    # Shorthand
    ctmech = mechanism.ctmech
    species_names = ctmech.species_names

    jacobian = np.zeros((len(species_names), len(species_names)))
    timescale = np.zeros((len(species_names)))

    ctmech.TPY = samples_dict['T'], samples_dict['P'], samples_dict['Y']

    concentration_save = np.zeros(len(ctmech.species()))
    concentration_new = np.zeros(len(ctmech.species()))
    old_source_term = np.zeros(len(ctmech.species()))
    new_source_term = np.zeros(len(ctmech.species()))

    concentration_save[:] = ctmech.concentrations[:]  # concentration kmol/m3
    old_source_term[:] = ctmech.net_production_rates[:]  # species source term kmol/m3/s

    if explicit:
        timescale = [abs(c / cdot) if c > abs_tol and cdot != 0 and abs(c / cdot) < 1 else 1.
                     for c, cdot in zip(concentration_save, old_source_term)]
    else:

        for index_spec, spec in enumerate(species_names):
            if concentration_save[index_spec] > abs_tol:
                dc = abs_tol + concentration_save[index_spec] * rel_tol  # perturbation in kmol/m3
                concentration_new[:] = concentration_save[:]
                concentration_new[index_spec] = dc + concentration_save[index_spec]  # only perturb one at a time
                ctmech.concentrations = concentration_new
                new_source_term[:] = ctmech.net_production_rates[:]
                jacobian[:, index_spec] = np.abs((new_source_term[:] - old_source_term[:])) / dc

                if concentration_weighting:
                    if jacobian[index_spec, index_spec] != 0 and jacobian[index_spec, index_spec] < 1e30:
                        timescale[index_spec] = concentration_save[index_spec] / jacobian[index_spec, index_spec]
                    else:
                        timescale[index_spec] = 1.0
                else:
                    if jacobian[index_spec, index_spec] != 0 and jacobian[index_spec, index_spec] < 1e30:
                        timescale[index_spec] = 1.0 / jacobian[index_spec, index_spec]
                    else:
                        timescale[index_spec] = 1.0

                if timescale[index_spec] > 1.0:
                    timescale[index_spec] = 1.0

            else:
                timescale[index_spec] = 1.0

    return timescale


def sort_species_by_solution(species, ctmech):
    """Sorts a list of species according to a Cantera Solution object order

    Parameters
    ----------
    species : list
        list of species to sort
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object

    Returns
    -------
    sorted_species = list
        list of species in the ctmech order


    """

    all_species = ctmech.species_names

    wrong_species = set(species) - set(all_species)

    if wrong_species:
        quit('Species ' + ', '.join(list(wrong_species)) + ' not in Cantera mechanism object')

    sorted_species = [all_species[all_species.index(spec)] for spec in all_species if spec in species]

    return sorted_species


def cutting_criterion(ranked_dict):
    """Computes criterion based on logarithmic mean

    Parameters
    ----------
    ranked_dict : dict
        directory on which applying the mean

    Returns
    -------
    criterion : `np.ndarray`
        criterion based on a logarithmic mean


    """

    # Cutting criterion estimation
    # Logarithmic mean of the values
    list_ranked_values = [ranked_dict[spec] for spec in ranked_dict if ranked_dict[spec]
                          if ranked_dict[spec] != 0]
    log_ranked_values = np.log(list_ranked_values)
    mean_log = np.mean(log_ranked_values)
    criterion = np.exp(mean_log)

    return criterion


def get_case_by_state(cases_list, reactor_type=None, P=-1, T=-1, phi=-1, fuel={}, oxidizer={}):
    """Get a list of cases corresponding to the specific initial state given

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    reactor_type : str
        type of reactor wanted (Default value = None)
    P : float
        initial pressure (Default value = -1)
    T : float
        initial temperature (Default value = -1)
    phi : float
        equivalence ratio :math:`\phi` (Default value = -1)
    fuel : dict
        fuel dictionary (Default value = {})
    oxidizer : dict
        oxidizer dictionary (Default value = {})

    Returns
    -------
    selected_cases : list
        list of class :func:`~ARCANE.cases.Case` objects corresponding to the given conditions


    """

    if reactor_type:
        cases_list_for_reactor = []

        found = False
        for reactor_name in kwdict.reactor_labels:
            for name in kwdict.reactor_labels[reactor_name]:
                if reactor_type == name:
                    reactor_type = reactor_name
                    found = True

        possible_names = [y for x in kwdict.reactor_labels.values() for y in x]

        if not found:
            logger.warning('The specified reactor type does not exist')
            logger.warning('Available reactor types are: ' + ', '.join(possible_names))

        types_of_reactors = []
        for case in cases_list:
            if reactor_type == case.reactor_type:
                cases_list_for_reactor.append(case)
            if case.reactor_type not in types_of_reactors:
                types_of_reactors.append(case.reactor_type)
    else:
        cases_list_for_reactor = cases_list
        types_of_reactors = [cases_list_for_reactor[0].reactor_type]

    number_of_input = 0
    if P > 0:
        number_of_input += 1
    if T > 0:
        number_of_input += 1
    if phi > 0:
        number_of_input += 1
    if fuel:
        number_of_input += 1
    if oxidizer:
        number_of_input += 1

    if number_of_input == 0 and len(types_of_reactors) == 1:
        logger.warning(
                'WARNING: the output will be the same as the list of cases as no initial state was specified ')
        logger.warning('You can specify P, T or phi with those keywords')

        selected_cases = cases_list

    else:

        selected_cases = []
        selected_ids = []

        for case in cases_list_for_reactor:
            pressure = case.pressure
            temperature = case.temperature
            if hasattr(case, 'phi'):
                phi_case = case.phi
            else:
                phi_case = None
            if hasattr(case, 'fuel'):
                fuel_case = case.fuel
            else:
                fuel_case = None
            if hasattr(case, 'oxidizer'):
                oxidizer_case = case.oxidizer
            else:
                oxidizer_case = None

            check = 0

            if pressure == P:
                check += 1
            if temperature == T:
                check += 1
            if phi_case == phi:
                check += 1
            if fuel_case == fuel:
                check += 1
            if oxidizer_case == oxidizer:
                check += 1

            if check >= number_of_input and case.myid not in selected_ids:
                selected_cases.append(case)
                selected_ids.append(case.myid)

    return selected_cases


def closest_case_id(case, mechanism, solution=False):
    """Checks in the cases directory what is the closest case to the current one

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism : :func:`~ARCANE.mechanisms.Mechanism`
        class :func:`~ARCANE.mechanisms.Mechanism` object
    solution : bool
        if True, check the files sol_***.xml and not the files profile_***

    Returns
    -------
    new_case_id : str
        case id that is the closest to the case (can be 'not_valid')


    """

    # Looking at the cases available in the cases directory
    computed_cases_id = database.get_cases_ids(mechanism, solution=solution)

    if database.database_system == 'database':
        computed_cases_id = [case_str for case_str in computed_cases_id if case_str.startswith(case.reactor_type)
                             and case_str != case.myid]
    else:
        computed_cases_id = [computed_case_id.replace(case.casedir + '/', '')
                             for computed_case_id in computed_cases_id
                             if computed_case_id.replace(case.casedir + '/', '') != case.myid
                             and computed_case_id.replace(case.casedir + '/', '').startswith(case.reactor_type)
                             and len(os.listdir(computed_case_id)) != 0]

    similarity_criterion = []

    for index, case_id in enumerate(computed_cases_id):
        meta_case = parse_id(case_id)

        if not meta_case:
            continue

        if meta_case.width != case.width or meta_case.soot_sections:
            similarity_criterion.append('0')
            continue

        criterion_values = []

        # Criterion on temperature
        if hasattr(case, 'temperature'):
            criterion = (min(case.temperature, meta_case.temperature) / max(case.temperature, meta_case.temperature)) **10
            criterion_values.append(round(criterion, 4))

        # Criterion on pressure (square root)
        if hasattr(case, 'pressure'):
            criterion = np.sqrt((min(case.pressure, meta_case.pressure))
                                / max(case.pressure, meta_case.pressure))
            criterion_values.append(round(criterion, 4))

        # Criterion on strain rate :math:`a` (square root)
        if hasattr(case, 'strain_rate'):
            criterion = np.sqrt((min(case.strain_rate, meta_case.strain_rate))
                                / max(case.strain_rate, meta_case.strain_rate))
            criterion_values.append(round(criterion, 4))

        # Criterion on scalar dissipation rate (square root)
        if hasattr(case, 'scalar_dissipation_rate'):
            criterion = np.sqrt((min(case.scalar_dissipation_rate, meta_case.scalar_dissipation_rate))
                                / max(case.scalar_dissipation_rate, meta_case.scalar_dissipation_rate))
            criterion_values.append(round(criterion, 4))

        # Criterion on equivalence ratio :math:`\phi` (cube)
        if hasattr(case, 'phi'):
            criterion = (min(case.phi, meta_case.phi)
                         / max(case.phi, meta_case.phi)) **3
            criterion_values.append(round(criterion, 4))

        # Criterion on fuel composition (square root)
        if hasattr(case, 'fuel') and isinstance(meta_case.fuel, dict):
            fuel_criterion_list = []
            for fuel_name in case.fuel:
                if fuel_name in meta_case.fuel:
                    fuel_criterion = np.sqrt(min(case.fuel[fuel_name], meta_case.fuel[fuel_name])
                                             / max(case.fuel[fuel_name], meta_case.fuel[fuel_name]))
                    fuel_criterion_list.append(fuel_criterion)
                else:
                    fuel_criterion = 0
                    fuel_criterion_list.append(fuel_criterion)
            criterion = np.sum(fuel_criterion_list)
            criterion_values.append(round(criterion / len(case.fuel), 4))

        # Criterion on oxidizer composition (square root)
        if hasattr(case, 'oxidizer') and isinstance(meta_case.fuel, dict):
            oxidizer_criterion_list = []
            for oxidizer_name in case.oxidizer:
                if oxidizer_name in meta_case.oxidizer:
                    oxidizer_criterion = np.sqrt(min(case.oxidizer[oxidizer_name], meta_case.oxidizer[oxidizer_name])
                                                 / max(case.oxidizer[oxidizer_name], meta_case.oxidizer[oxidizer_name]))
                    oxidizer_criterion_list.append(oxidizer_criterion)
                else:
                    oxidizer_criterion = 0
                    oxidizer_criterion_list.append(oxidizer_criterion)
            criterion = np.sum(oxidizer_criterion_list)
            criterion_values.append(round(criterion / len(case.oxidizer), 4))

        # Criterion on composition (square root)
        if meta_case.composition:
            composition_criterion_list = []
            for composition_name in case.composition:
                if composition_name in meta_case.composition:
                    composition_criterion = np.sqrt(
                            min(case.composition[composition_name], meta_case.composition[composition_name])
                            / max(case.composition[composition_name], meta_case.composition[composition_name]))
                    composition_criterion_list.append(composition_criterion)
                else:
                    composition_criterion = 0
                    composition_criterion_list.append(composition_criterion)
            criterion = np.sum(composition_criterion_list)
            criterion_values.append(round(criterion / len(case.composition), 4))

        # WTF type bug changing the value of the floats
        # the following lines are dirty but working
        similarity_criterion.append(str(np.prod(criterion_values)))

    similarity_criterion = [float(value) for value in similarity_criterion]

    # Selecting the case with the highest value
    new_case_id = 'not_valid'

    if similarity_criterion:
        #
        # zipped_lists = zip(similarity_criterion, computed_cases_id)
        # sorted_pairs = sorted(zipped_lists)
        #
        # tuples = zip(*sorted_pairs)
        # list1, list2 = [list(tuple) for tuple in tuples]
        #
        # for index, crit in enumerate(list1):
        #     print(crit, list2[index])

        best_criterion = np.max(similarity_criterion)
        if best_criterion > 0:
            index_case = similarity_criterion.index(best_criterion)
            logger.debug('Found a case with ' + str(best_criterion) + ' similarity')
            if best_criterion > 0.5:
                new_case_id = computed_cases_id[index_case]
            else:
                logger.debug('Not good enough')
        else:
            logger.debug('No similar case found')

    return new_case_id


def case_from_id(id_string, mechanism):
    """Creates a case from the an id

    Parameters
    ----------
    id_string : str
        id of the case
    mechanism : :func:`~ARCANE.mechanisms.Mechanism`
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object corresponding to the given id


    """

    meta_case = parse_id(id_string)

    case = cases.create_case(reactor=meta_case.reactor_type,
                             mechanism=mechanism,
                             fuel=meta_case.fuel,
                             oxidizer=meta_case.oxidizer,
                             pressure=meta_case.pressure,
                             temperature=meta_case.temperature,
                             phi=meta_case.phi)[0]

    return case


def parse_id(id_string):
    """Get the basics parameters of the case from the id

    Parameters
    ----------
    id_string : str
        id of the case

    Returns
    -------
    meta_case : class :func:`~ARCANE.tools.parse_id.MetaCase` object
        class :func:`~ARCANE.tools.parse_id.MetaCase` object containing the case basics parameters


    """

    class MetaCase(object):
        """MetaCase class, containing all class-nonspecific methods and attributes"""

        # Generator
        def __init__(self):
            """Constructor of class MetaCase"""
            self.myid = None
            self.temperature = None
            self.fuel_temperature = None
            self.oxidizer_temperature = None
            self.pressure = None
            self.phi = 0
            self.composition = {}
            self.reactor_type = None
            self.fuel = {}
            self.oxidizer = {}
            self.strain_rate = None
            self.scalar_dissipation_rate = None
            self.soot_sections = 0
            self.soot_precursors = None
            self.width = None

    meta_case = MetaCase()
    meta_case.myid = id_string

    # Delimiters
    delimiters_list = ['p', 'Tfuel', 'Toxi', 'T',
                       'Compo', 'phi',
                       'Fuel', 'Oxi',
                       'a', 'chist', '$',
                       'SOOT', 'PREC', 'width']

    cuts = {}
    trunc = meta_case.myid

    used_delimiters = []

    try:
        count = 1
        for index, delimiter in enumerate(delimiters_list):
            if delimiter in trunc:
                cut = trunc.split(delimiter, 1)
                if delimiter == delimiters_list[0]:
                    cuts['reactor_type'] = cut[0]
                    trunc = cut[1]
                else:
                    cuts[delimiters_list[index - count]] = cut[0]
                    trunc = cut[1]
                last_delimiter = delimiter
                used_delimiters.append(delimiter)
                count = 1
            else:
                count += 1

        cuts[last_delimiter] = trunc

        if 'Fuel' in used_delimiters:
            fuel_list = cuts['Fuel'].split('%')
            index = 1
            while index < len(fuel_list) - 1:
                if len(fuel_list[index]) > 3:
                    new_split = [fuel_list[index][:3], fuel_list[index][3:]]
                    fuel_list.remove(fuel_list[index])
                    fuel_list.insert(index, new_split[0])
                    fuel_list.insert(index + 1, new_split[1])
                    index += 2

            fuel = {}
            for index in range(0, len(fuel_list), 2):
                fuel[fuel_list[index]] = float(fuel_list[index + 1]) / 100

            cuts['Fuel'] = fuel

        if 'Oxi' in used_delimiters:
            oxidizer_list = cuts['Oxi'].split('%')
            index = 1
            while index < len(oxidizer_list) - 1:
                if len(oxidizer_list[index]) > 3:
                    new_split = [oxidizer_list[index][:3], oxidizer_list[index][3:]]
                    oxidizer_list.remove(oxidizer_list[index])
                    oxidizer_list.insert(index, new_split[0])
                    oxidizer_list.insert(index + 1, new_split[1])
                    index += 2

            oxidizer = {}
            for index in range(0, len(oxidizer_list), 2):
                oxidizer[oxidizer_list[index]] = float(oxidizer_list[index + 1]) / 100

            cuts['Oxi'] = oxidizer

        if 'Compo' in used_delimiters:
            composition_list = cuts['Compo'].split('%')
            index = 1
            while index < len(composition_list) - 1:
                if len(composition_list[index]) > 3:
                    new_split = [composition_list[index][:3], composition_list[index][3:]]
                    composition_list.remove(composition_list[index])
                    composition_list.insert(index, new_split[0])
                    composition_list.insert(index + 1, new_split[1])
                    index += 2

            composition = {}
            for index in range(0, len(composition_list), 2):
                composition[composition_list[index]] = float(composition_list[index + 1]) / 100

            cuts['Compo'] = composition
        
        if 'width' not in used_delimiters and cuts['reactor_type'] in kwdict.reactor_default_values:
            cuts['width'] = kwdict.reactor_default_values[cuts['reactor_type']]['width'] * 100

        # Setting attributes of the meta case
        attribute_cuts_correspondence = {}
        attribute_cuts_correspondence['reactor_type'] = 'reactor_type'
        attribute_cuts_correspondence['p'] = 'pressure'
        attribute_cuts_correspondence['T'] = 'temperature'
        attribute_cuts_correspondence['Tfuel'] = 'fuel_temperature'
        attribute_cuts_correspondence['Toxi'] = 'oxidizer_temperature'
        attribute_cuts_correspondence['phi'] = 'phi'
        attribute_cuts_correspondence['Fuel'] = 'fuel'
        attribute_cuts_correspondence['Oxi'] = 'oxidizer'
        attribute_cuts_correspondence['Compo'] = 'composition'
        attribute_cuts_correspondence['a'] = 'strain_rate'
        attribute_cuts_correspondence['chist'] = 'scalar_dissipation_rate'
        attribute_cuts_correspondence['width'] = 'width'
        attribute_cuts_correspondence['SOOT'] = 'soot_sections'

        # Converting to floats
        for delimiter in cuts:
            if delimiter in ['T', 'Tfuel', 'Toxi', 'a', 'chist', 'SOOT']:
                cuts[delimiter] = float(cuts[delimiter])
            elif delimiter in ['phi', 'width']:
                cuts[delimiter] = float(cuts[delimiter]) / 100
            elif delimiter in ['p']:
                cuts[delimiter] = float(cuts[delimiter]) * 1000

            setattr(meta_case, attribute_cuts_correspondence[delimiter], cuts[delimiter])

    except ValueError:
        logger.debug('Non-standard naming of the case; cannot be parsed')
        meta_case = None

    return meta_case


def clean(dir='all'):
    """Cleans the directories in which data has been written

    Parameters
    ----------
    dir : str
        name of the directory to clean (Default value = 'all')

    """
    if dir == 'all':
        shutil.rmtree('cases')
        shutil.rmtree('mechs')
        shutil.rmtree('reduced')
    elif dir == 'cases':
        shutil.rmtree('cases')
    elif dir == 'mechs':
        shutil.rmtree('mechs')
    elif dir == 'reduced':
        shutil.rmtree('reduced')
    else:
        shutil.rmtree(dir)


def format_to_significant(x, n):
    """Rounds a number to 15 significant digits to be coherent with fortran precision

    Parameters
    ----------
    x :
        number to round
    n :
        number of significant digits

    Returns
    -------
     : str
        rounded number at 15 significant digits


    """

    # Python number display overcome

    # return format(x, '.' + str(n) + 'e')
    return "{:e}".format(x)


def compute_Schmidt_numbers(ctmech, spec_list='all'):
    """Computing Schmidt number for the input species

    Parameters
    ----------
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object
    spec_list : list
        list of species to retrieve (default is all species)

    Returns
    -------
    Sch : list
        list of Schmidt corresponding to the given species


    """

    if ctmech.transport_model == 'Transport':
        ctmech.transport_model = 'Mix'

    if spec_list == 'all':
        spec_list = ctmech.species_names

    Sch = []
    for index, spec in enumerate(spec_list):
        index_spec = ctmech.species_index(spec)
        if (ctmech.mix_diff_coeffs[index_spec] * ctmech.density) != 0:
            Sch.append(ctmech.viscosity / (ctmech.mix_diff_coeffs[index_spec] * ctmech.density))
        else:
            Sch.append(np.nan)

    return Sch


def Lewis_profile(case, mechanism, spec):
    """Computes a species Lewis number profile on a 1D case

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object
    spec : str
        species name


    Returns
    -------
    Le_profile : list
        profile of Lewis along the case grid


    """
    data_dict = case.data_dict(mechanism)
    grid = data_dict['Grid']

    ctmech = mechanism.ctmech
    ctmech.transport_model = 'Mix'

    spec_index = ctmech.species_index(spec)

    Le_profile = []
    for index, value in enumerate(grid):
        ctmech.TPY = data_dict['T'][index], data_dict['P'][index], data_dict['Y'][index]
        Le_profile.append(compute_Schmidt_numbers(ctmech)[spec_index] / compute_Prandtl_number(ctmech))

    return Le_profile


def Schmidt_profile(case, mechanism, spec):
    """Computes a species Schmidt number profile on a 1D case

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object
    spec : str
        species name


    Returns
    -------
    Sc_profile : list
        profile of Schmidt along the case grid


    """
    data_dict = case.data_dict(mechanism)
    grid = data_dict['Grid']

    ctmech = mechanism.ctmech
    ctmech.transport_model = 'Mix'

    spec_index = ctmech.species_index(spec)

    Sc_profile = []
    for index, value in enumerate(grid):
        ctmech.TPY = data_dict['T'][index], data_dict['P'][index], data_dict['Y'][index]
        Sc_profile.append(compute_Schmidt_numbers(ctmech)[spec_index])

    return Sc_profile


def compute_Prandtl_number(ctmech):
    """Compute mixture Prandtl number

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object


    Returns
    -------
    Pr : float
        Prandtl Number


    """

    if ctmech.transport_model == 'Transport':
        ctmech.transport_model = 'Mix'

    Pr = ctmech.viscosity * ctmech.cp_mass / ctmech.thermal_conductivity

    return Pr


def Prandtl_profile(case, mechanism):
    """Computes a Prandtl number profile on a 1D case

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object
    spec : str
        species name

    Returns
    -------
    Pr_profile : list
        profile of Prandtl along the case grid


    """
    data_dict = case.data_dict(mechanism)
    grid = data_dict['Grid']

    ctmech = mechanism.ctmech
    ctmech.transport_model = 'Mix'

    Pr_profile = []
    for index, value in enumerate(grid):
        ctmech.TPY = data_dict['T'][index], data_dict['P'][index], data_dict['Y'][index]
        Pr_profile.append(compute_Prandtl_number(ctmech))

    return Pr_profile


def get_qss_species(cti):
    """Get the list of QSS species from a formatted cti file

    Parameters
    ----------
    cti : str
        Cantera input file

    Returns
    -------
    species_qss_names : list
        list of QSS species


    """

    file = open(cti)
    lines = file.readlines()

    of_interest = ''

    for index, line in enumerate(lines):
        if 'species_qss' in line:
            index_quote = line.find('"""')
            of_interest = line[index_quote + 3:]
            break

    while of_interest:
        index += 1
        if lines[index].startswith('#'):
            of_interest += lines[index]
        else:
            break

    of_interest = of_interest.replace('""", ', '')
    go_on = True
    ref = of_interest
    while go_on:
        output = of_interest.replace('  ', '')
        if output != ref:
            ref = output
        else:
            go_on = False
    of_interest = of_interest.replace('# ', ' ')
    of_interest = of_interest.replace('#', ' ')
    of_interest = of_interest.replace(' \n', '')
    of_interest = of_interest.replace('\n', '')
    of_interest = of_interest.replace(' ', "', '")
    of_interest = of_interest.replace('""",', "")
    if not of_interest.endswith("', '"):
        of_interest += "', '"
    of_interest = '[' + of_interest[3:-3] + ']'
    if of_interest:
        species_qss_names = eval(of_interest)
    else:
        species_qss_names = []

    species_qss_names = [spec for spec in species_qss_names if spec != '']

    file.close()

    return species_qss_names


def get_kinetics(cti):
    """Extracts the type of kinetics
    (very simplistic for now)

    Parameters
    ----------
    cti : str
        Cantera input file

    Returns
    -------
    kinetics : str
        string stating the type of kinetics ('custom' or None)

    """

    file = open(cti)
    lines = file.readlines()

    kinetics = None

    for index, line in enumerate(lines):
        transformed_line = line.replace(' ', '').replace("'", '').replace('"', '').replace(',', '').replace('\n', '')
        if transformed_line == 'kinetics=custom':
            kinetics = 'custom'
            break

        elif all([word in transformed_line for word in ['pecies', 'ata']]) \
                or all([word in transformed_line for word in ['name', 'species']]):
            break

    return kinetics


def get_number_of_reactions(cti):
    """Extracts the number of reactions from a formatted cti file

    Parameters
    ----------
    cti : str
        Cantera input file

    Returns
    -------
    nr : int
        Number of reactions


    """
    for line in reversed(open(cti).readlines()):
        if 'Reaction' in line:
            number_string = line.replace('##  Reaction ', '')
            break

    nr = int(number_string)

    return nr


def skeletal_from_custom_cti(cti, file_name='skeletal.cti'):
    """Retrieve the cti file for the skeletal mechanism from a reduced one

    Parameters
    ----------
    cti : str
        Cantera input file
    file_name : str
        output cti name (Default value = 'skeletal.cti')


    """

    f = open(file_name, 'w')

    cti_file = open(cti)
    file = cti_file.readlines()

    spec_part = True
    qss_part = False
    reactions_part = False
    skip_counter = 0

    for index, line in enumerate(file):
        if skip_counter > 0:
            skip_counter -= 1
        else:
            if '""",' in line and spec_part:
                new_line = line.replace('""",', '')
                spec_part = False
                qss_part = True
                f.write(new_line)

            elif 'species_qss' in line:

                line = line.replace('""",', "Thiswillbereplaced")
                split_line = line.split('"""')
                first_part = split_line[0]
                second_part = split_line[1]
                new_line = first_part.replace(' ', '').replace('  ', '').replace('#', '                   ')
                new_line = new_line.replace('species_qss', '').replace('=', '').replace('"""', '')
                new_line = new_line + second_part
                new_line = new_line.replace("Thiswillbereplaced", '""",')

                f.write(new_line)

            elif line.startswith('#') and (qss_part or reactions_part):
                new_line = line.replace('##', '?')
                new_line = new_line.replace('#', '')
                new_line = new_line.replace('?', '#')
                f.write(new_line)

            elif line.replace(' ', '').startswith('reactions='):
                qss_part = False
                f.write(line)

            elif line.replace(' ', '').replace("'", '"').startswith('kinetics="custom",'):
                continue

            elif line.replace(' ', '').replace("'", '"').startswith('transport="AVBP",'):
                 continue

            elif line.replace(' ', '').startswith('#Dummyreaction'):
                skip_counter = 2
                reactions_part = True
                continue
            else:
                f.write(line)

    f.close()
    cti_file.close()


def copy_reaction(reaction):
    """Copying manually the reactions as they are not pickable

    Parameters
    ----------
    reaction : `Cantera.Reaction`
        class `Cantera.Reaction` object


    Returns
    -------
    new_reaction : class `Cantera.Reaction` object
        new class `Cantera.Reaction` object, can be either `Cantera.ElementaryReaction`, `Cantera.ThreeBodyReaction`,
        `Cantera.FalloffReaction`, `Cantera.PlogReaction`


    """

    # Extracting data from the old reaction
    reactants = reaction.reactants
    products = reaction.products

    reaction_type = reaction.reaction_type

    if reaction_type == 1:

        # New elementary reaction
        new_reaction = ct.ElementaryReaction(reactants, products)
        new_reaction.rate = reaction.rate

    elif reaction_type == 2:

        # New three bodies reaction
        new_reaction = ct.ThreeBodyReaction(reactants, products)
        new_reaction.rate = reaction.rate

        new_reaction.efficiencies = reaction.efficiencies

    elif reaction_type == 4:

        # New falloff reaction
        new_reaction = ct.FalloffReaction(reactants, products)

        new_reaction.high_rate = reaction.high_rate
        new_reaction.low_rate = reaction.low_rate

        new_reaction.efficiencies = reaction.efficiencies

        new_reaction.falloff = reaction.falloff

    elif reaction_type == 5:

        # New pressure dependent reaction
        new_reaction = ct.PlogReaction(reactants, products)

        new_reaction.rates = reaction.rates

    else:
        logger.error('The copy of this type of reaction is not implemented thus the lumping cannot be done.')

    new_reaction.ID = reaction.ID
    new_reaction.duplicate = reaction.duplicate
    new_reaction.orders = reaction.orders
    new_reaction.reversible = reaction.reversible

    return new_reaction


def dual_flow_equivalence_ratio(ctmech, first_mixture_dict, second_mixture_dict):
    """Computes the global equivalence ratio :math:`\phi_g` for a configuration with multiple fuels in diffusion case

    Parameters
    ----------
    ctmech : `Cantera.Solution`
        class `Cantera.Solution` object
    first_mixture_dict : dict
        Dictionary with mixture infos in the following format:

        ``first_mixture_dict = {'fuel': 'CH4',``

        ``'oxidizer': {'O2': 0.21, 'N2': 0.79},``

        ``'phi': 0.8, 'velocity': 2}``
    second_mixture_dict : dict
        same as ``first_mixture_dict`` for the ``second_mixture_dict``

    Returns
    -------
     : float
        Global equivalence ratio :math:`\phi_g`


    """

    # Checking that every parameter is there
    # Checking first dictionary
    fuel_first = False
    oxidizer_first = False
    phi_first = False
    velocity_first = False

    if 'fuel' in first_mixture_dict:
        fuel_first = True
    if 'oxidizer' in first_mixture_dict:
        oxidizer_first = True
    if 'phi' in first_mixture_dict:
        phi_first = True
    if 'velocity' in first_mixture_dict:
        velocity_first = True

    # Checking second dictionary
    fuel_second = False
    oxidizer_second = False
    phi_second = False
    velocity_second = False

    if 'fuel' in second_mixture_dict:
        fuel_second = True
    if 'oxidizer' in second_mixture_dict:
        oxidizer_second = True
    if 'phi' in second_mixture_dict:
        phi_second = True
    if 'velocity' in second_mixture_dict:
        velocity_second = True

    # Raising errors if parameters are missing
    if phi_first and not fuel_first and not oxidizer_first:
        logger.error('You cannot specify an equivalence ratio :math:`\phi` without specifying both fuel and oxidizer')
        logger.error('--->' + str(first_mixture_dict))
        quit()
    if not phi_first and fuel_first and oxidizer_first:
        logger.error(
                "You need to specify an equivalence ratio :math:`\phi` for the first mixture with the keyword 'phi' if fuel and oxidizer are both specified")
        quit()

    if phi_second and not fuel_second and not oxidizer_second:
        logger.error('You cannot specify an equivalence ratio :math:`\phi` without specifying both fuel and oxidizer')
        logger.error('--->' + str(second_mixture_dict))
        quit()
    if not phi_second and fuel_second and oxidizer_second:
        logger.error(
                "You need to specify an equivalence ratio :math:`\phi` for the second mixture with the keyword 'phi' if fuel and oxidizer are both specified")
        quit()

    if not velocity_first:
        logger.error("Please specify a velocity for the first mixture with the keyword 'velocity'")
        quit()

    if not velocity_second:
        logger.error("Please specify a velocity for the second mixture with the keyword 'velocity'")
        quit()

    # Converting strings to directory
    # And extracting composition of the molecules from Cantera mechanism
    composition_dict = {}
    for which_dict in [first_mixture_dict, second_mixture_dict]:
        for mix in ['fuel', 'oxidizer']:
            if mix in which_dict:
                if type(which_dict[mix]) == str:
                    which_dict[mix] = {which_dict[mix]: 1}
                # Values normalization
                min_value = min(which_dict[mix].values())
                for spec in which_dict[mix]:
                    which_dict[mix][spec] /= min_value
                    if spec not in composition_dict:
                        composition_dict[spec] = {}
                        for element in ctmech.element_names:
                            composition_dict[spec][element] = ctmech.n_atoms(spec, element)

    # Velocity ratio
    velocity_ratio = first_mixture_dict['velocity'] / second_mixture_dict['velocity']

    # Equilibrate equation for first mixture
    if phi_first:
        # Compute phi for oxidizer of type CxHyOz
        soxi_num = 0.0
        for key in first_mixture_dict['fuel']:
            soxi_num += first_mixture_dict['fuel'][key] * (- 4 * composition_dict[key]['C']
                                                           - composition_dict[key]['H']
                                                           + 2 * composition_dict[key]['O'])

        # Compute phi
        soxi_denom = 0.0
        for key in first_mixture_dict['oxidizer']:
            soxi_denom += first_mixture_dict['oxidizer'][key] * (4 * composition_dict[key]['C']
                                                                 + composition_dict[key]['H']
                                                                 - 2 * composition_dict[key]['O'])

        soxi_first = soxi_num / soxi_denom

    elif fuel_first:
        soxi_first = 0

    elif oxidizer_first:
        soxi_first = 1

    # Equilibrate equation for first mixture
    if phi_second:
        # Compute phi for oxidizer of type CxHyOz
        soxi_num = 0.0
        for key in second_mixture_dict['fuel']:
            soxi_num += second_mixture_dict['fuel'][key] * (- 4 * composition_dict[key]['C']
                                                            - composition_dict[key]['H']
                                                            + 2 * composition_dict[key]['O'])

        # Compute phi
        soxi_denom = 0.0
        for key in second_mixture_dict['oxidizer']:
            soxi_denom += second_mixture_dict['oxidizer'][key] * (4 * composition_dict[key]['C']
                                                                  + composition_dict[key]['H']
                                                                  - 2 * composition_dict[key]['O'])

        soxi_second = soxi_num / soxi_denom

    elif fuel_second:
        soxi_second = 0

    elif oxidizer_second:
        soxi_second = 1

    # Global oxidizer stoechiometric coefficient
    global_soxi = velocity_ratio * soxi_first + soxi_second

    new_first = first_mixture_dict.copy()
    new_second = second_mixture_dict.copy()

    full_mixture = {'fuel': {}, 'oxidizer': {}}
    # Retrieve equivalent fuel
    if 'fuel' in new_first:
        for spec in new_first['fuel']:
            new_first['fuel'][spec] *= velocity_ratio
            if phi_first:
                new_first['fuel'][spec] *= new_first['phi']
            if spec not in full_mixture['fuel']:
                full_mixture['fuel'][spec] = new_first['fuel'][spec]
            else:
                full_mixture['fuel'][spec] += new_first['fuel'][spec]

    if 'fuel' in new_second:
        for spec in new_second['fuel']:
            if phi_second:
                new_second['fuel'][spec] *= new_second['phi']
            if spec not in full_mixture['fuel']:
                full_mixture['fuel'][spec] = new_second['fuel'][spec]
            else:
                full_mixture['fuel'][spec] += new_second['fuel'][spec]

    # Retrieve equivalent oxidizer
    if 'oxidizer' in new_first:
        for spec in new_first['oxidizer']:
            new_first['oxidizer'][spec] *= velocity_ratio
            if spec not in full_mixture['oxidizer']:
                full_mixture['oxidizer'][spec] = new_first['oxidizer'][spec]
            else:
                full_mixture['oxidizer'][spec] += new_first['oxidizer'][spec]

    if 'oxidizer' in new_second:
        for spec in new_second['oxidizer']:
            if spec not in full_mixture['oxidizer']:
                full_mixture['oxidizer'][spec] = new_second['oxidizer'][spec]
            else:
                full_mixture['oxidizer'][spec] += new_second['oxidizer'][spec]

    # Normalizing
    min_value = min(full_mixture['oxidizer'].values())
    for spec in full_mixture['oxidizer']:
        full_mixture['oxidizer'][spec] /= min_value

    # Equilibrate the global equation
    soxi_num = 0.0
    for key in full_mixture['fuel']:
        soxi_num += full_mixture['fuel'][key] * (- 4 * composition_dict[key]['C']
                                                 - composition_dict[key]['H']
                                                 + 2 * composition_dict[key]['O'])

    soxi_denom = 0.0
    for key in full_mixture['oxidizer']:
        soxi_denom += full_mixture['oxidizer'][key] * (4 * composition_dict[key]['C']
                                                       + composition_dict[key]['H']
                                                       - 2 * composition_dict[key]['O'])

    global_soxi_st = soxi_num / soxi_denom

    global_equivalence_ratio = global_soxi_st / global_soxi

    return round(global_equivalence_ratio, 4)


def Arrhenius_fit(temp_list, arrhenius_values, option='full'):
    """Fitting Arrhenius parameters on a given dataset

    Parameters
    ----------
    temp_list :
        list of temperature
    arrhenius_values :
        list of reaction coefficient values
    option :
        specifies which type of function will be used (Arrhenius, Arrhenius without T dependence or constant) (Default value = 'full')

    Returns
    -------
    new_coefficients : list
        list of coefficient for the Arrhenius


    """
    import warnings
    warnings.filterwarnings('ignore')
    A = 1.0
    beta = 0.0
    Ea = 0.0
    params = [A, beta, Ea]

    temperature = np.array(temp_list)
    # ln_arrhenius = np.array([np.log(Arrhenius_values[i]) for i in range(len(temp_list))])
    k_values = np.array([arrhenius_values[i] for i in range(len(temp_list))])

    error_list = []
    list_of_coefficients = []

    if option == 'best':
        options_to_try = ['full', 'simplified', 'constant']
    else:
        options_to_try = [option]

    for option_local in options_to_try:
        new_coefficients = []
        arguments = (temperature, k_values, option_local)

        # LEASTSQ FIT
        display.display(False)
        fit = optimize.leastsq(Arrhenius_leastsq_function, params, args=arguments)
        display.display(True)

        for j in range(3):
            new_coefficients.append(fit[0][j])

        list_of_coefficients.append(new_coefficients)
        # Testing error
        error = Arrhenius_leastsq_function(new_coefficients, *arguments)
        error_list.append(np.sqrt(sum([x ** 2 for x in error])))

    index_min = error_list.index(min(error_list))

    new_coefficients = list_of_coefficients[index_min]

    if option == 'best':
        logger.debug('Arrhenius has been fitted using ' + options_to_try[index_min] + ' function\n')

    return new_coefficients


def Arrhenius_leastsq_function(params, *args):
    """Least square function on Arrhenius function

    Parameters
    ----------
    params :
        list of coefficients to initialise the fit (A, beta, Ea)
    *args :
        parameters to fit (temperature, values of Arrhenius and option for specifying the function to use)


    Returns
    -------
     : ndarray
        value to minimize for the fitting of the Arrhenius

    """

    A = params[0]
    beta = params[1]
    Ea = params[2]

    T = args[0]
    k = args[1]
    option = args[2]

    kfit = np.empty(T.shape)

    kfit[:] = Arrhenius_function(T, A, beta, Ea, option=option)

    return k - kfit


def Arrhenius_function(T, A, beta, Ea, option='full'):
    """Computing the Arrhenius function (modified or not according to option)

    Parameters
    ----------
    T :
        Temperature
    A :
        pre-exponential coefficient
    beta :
        temperature exponent
    Ea :
        activation energy
    option :
        specifies which type of function will be used (Arrhenius, Arrhenius with T dependence or constant) (Default value = 'full')

    Returns
    -------
    k : float
        Arrhenius coefficient


    """

    R = 8314.4621

    if option == 'constant':
        k = A
    elif option == 'simplified':
        k = A * np.exp(- Ea / (R * T))
    else:
        k = A * T ** beta * np.exp(- Ea / (R * T))

    return k


def calculate_lnk(T, ln_A, beta, Ea_R):
    """Formula of ln(k) using Arrhenius expression, use SI units

    Parameters
    ----------
    T :
        temperature
    ln_A :
        natural log of pre-exponential factor
    beta :
        temperature exponent
    Ea_R :
        activation energy divided by ideal gas constant

    Returns
    -------
     : float
        value of ln(k)


    """
    return ln_A + (beta * np.log(T)) + (-Ea_R / T)


def calculate_lnk_simple(T, ln_A, Ea_R):
    """Formula of ln(k) using Arrhenius expression, use SI units

    Parameters
    ----------
    T :
        temperature
    ln_A :
        natural log of pre-exponential factor
    Ea_R :
        activation energy divided by ideal gas constant

    Returns
    -------
     : float
        value if ln(k) without beta taken into account



    """
    return ln_A + (-Ea_R / T)


def change_rate(reaction, A_factor=1.0, A_addition=0.0,
                b_factor=1.0, b_addition=0.0,
                Ea_factor=1.0, Ea_addition=0.0):
    """Function multiplying the pre-exponential factor of a reaction by a given factor
    or adding a to it

    Parameters
    ----------
    reaction :
        class `Cantera.Reaction` object to modify
    A_factor :
        factor by which the pre-exponential factor will be multiplied (Default value = 1.0)
    A_addition :
        value added to the pre-exponential factor (Default value = 0.0)
    b_factor :
        factor by which the temperature exponent will be multiplied (Default value = 1.0)
    b_addition :
        value added to the temperature exponent (Default value = 0.0)
    Ea_factor :
        factor by which the activation energy will be multiplied (Default value = 1.0)
    Ea_addition :
        value added to the activation energy (Default value = 0.0)

    Returns
    -------
    reaction : class `Cantera.Reaction` object
        modified class `Cantera.Reaction` object


    """

    if reaction.reaction_type in [1, 2]:

        rate = reaction.rate
        A = rate.pre_exponential_factor
        b = rate.temperature_exponent
        Ea = rate.activation_energy

        reaction.rate = ct.Arrhenius(A * A_factor + A_addition,
                                     b * b_factor + b_addition,
                                     Ea * Ea_factor + Ea_addition)

    elif reaction.reaction_type == 4:

        high_rate = reaction.high_rate
        A = high_rate.pre_exponential_factor
        b = high_rate.temperature_exponent
        Ea = high_rate.activation_energy

        reaction.high_rate = ct.Arrhenius(A * A_factor + A_addition,
                                          b * b_factor + b_addition,
                                          Ea * Ea_factor + Ea_addition)

        low_rate = reaction.low_rate
        A = low_rate.pre_exponential_factor
        b = low_rate.temperature_exponent
        Ea = low_rate.activation_energy

        reaction.low_rate = ct.Arrhenius(A * A_factor + A_addition,
                                         b * b_factor + b_addition,
                                         Ea * Ea_factor + Ea_addition)

    return reaction


def change_reaction_by_species(ctmech, species_names_list, A=1, b=0, Ea=0, new_names_list=''):
    """Modifies all the reactions that contains the specified species

    Parameters
    ----------
    ctmech :
        class `Cantera.Solution` object
    species_names_list :
        name or list of names of the targeted species
    A :
        pre-exponential factor (Default value = 1)
    b :
        temperature exponent (Default value = 0)
    Ea :
        activation energy (Default value = 0)
    new_names_list :
        new name or list of new names of the species (if to be changed) (Default value = '')

    Returns
    -------
    new_reactions : list
        list of all the class `Cantera.Reaction` objects corresponding to the given species


    """

    species_names = ctmech.species_names
    reactions = [copy_reaction(reac) for reac in ctmech.reactions()]
    new_reactions = reactions

    nur = ct.Kinetics.reactant_stoich_coeffs(ctmech)
    nup = ct.Kinetics.product_stoich_coeffs(ctmech)

    if not isinstance(species_names_list, list):
        species_names_list = [species_names_list]

    if not isinstance(new_names_list, list):
        new_names_list = [new_names_list]

    if new_names_list:
        if not len(species_names_list) == len(new_names_list):
            logger.error('The list of new names must be the same length as the list of species')

    for index, species_name in enumerate(species_names_list):

        if species_name in species_names:

            # Global index of the species to lump
            global_index = species_names.index(species_name)

            new_name = new_names_list[index]

            reactions_to_discard = list()

            # Indices of reactions containing the species
            spec_reac_indices = list(np.nonzero(nur[global_index, :])[0])
            spec_prod_indices = list(np.nonzero(nup[global_index, :])[0])

            for reac_index in spec_reac_indices:

                # Get reactants dictionary
                reactants_dict = reactions[reac_index].reactants

                # Get stoechimoetric coefficient
                stoech_spec = reactants_dict[species_name]

                if new_name:
                    # Checking if new species already in reaction,
                    # if that's case the number of molecules are added
                    if new_name in reactants_dict:
                        stoech_new_spec = reactants_dict[new_name]
                    else:
                        stoech_new_spec = 0

                    # Remove previous species
                    reactants_dict.pop(species_name)

                    reactants_dict[new_name] = stoech_new_spec + stoech_spec

                # Set reactants dictionary
                reactions[reac_index].reactants = reactants_dict

                # Changing the rate of the reaction
                reactions[reac_index] = change_rate(reactions[reac_index],
                                                    A_factor=A, b_addition=b, Ea_addition=Ea)

            for reac_index in spec_prod_indices:

                # Get reactants dictionary
                products_dict = reactions[reac_index].products

                # Get stoechimoetric coefficient
                stoech_spec = products_dict[species_name]

                if new_name:
                    # Checking if new species already in reaction,
                    # if that's case the number of molecules are added
                    if new_name in products_dict:
                        stoech_new_spec = products_dict[new_name]
                    else:
                        stoech_new_spec = 0
                    # Remove previous species
                    products_dict.pop(species_name)

                    # Add lumped one
                    products_dict[new_name] = stoech_new_spec + stoech_spec

                # Set reactants dictionary
                reactions[reac_index].products = products_dict

            # Discard isomerization reactions that are becoming A <=> A
            for reac_index, reac in enumerate(reactions):
                if reac.products == reac.reactants and ctmech.is_reversible(reac_index):
                    reactions_to_discard.append(reac_index)

                if reac.reaction_type in [2, 4]:
                    new_efficiencies = reac.efficiencies
                    if species_name in reac.efficiencies:
                        new_efficiencies.pop(species_name)
                    reac.efficiencies = new_efficiencies

            # Retrieve all reactions except redundant ones
            new_reactions = [reac for index_reac, reac in enumerate(reactions)
                             if index_reac not in reactions_to_discard]

            reactions = [copy_reaction(reac) for reac in new_reactions]

    return new_reactions


def get_thermo(mechanism, temperature, composition_dict={}, option=''):
    """Function to get all the thermal variables for a mixture

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object with the species to lump
    temperature :
        An array of Temperature on which we have to estimate all the thermal variables (heat capacity cp, enthalpy H and entropy S)
    composition_dict :
        Composition of the gas to estimate all the the thermal variables of the mixture.
        If not set, it takes the default value of the mechanism
    option :
        (Default value = '')

    """

    ctmech = mechanism.ctmech
    if composition_dict:
        mechanism.ctmech.X = composition_dict

    pressure = ctmech.species()[0].thermo.reference_pressure

    mechanism.ctmech.TP = temperature, pressure
    if option == 'cp':
        cp = ctmech.cp * ctmech.mean_molecular_weight
        return cp

    elif option == 'enthalpy':
        enthalpy = ctmech.HP[0] * ctmech.mean_molecular_weight
        return enthalpy

    elif option == 'entropy':
        entropy = ctmech.SP[0] * ctmech.mean_molecular_weight
        return entropy

    elif option == 'all':
        cp = ctmech.cp * ctmech.mean_molecular_weight
        enthalpy = ctmech.HP[0] * ctmech.mean_molecular_weight
        entropy = ctmech.SP[0] * ctmech.mean_molecular_weight
        return cp, enthalpy, entropy


def NASA(temperature, a1, a2, a3, a4, a5, a6, a7, quantity):
    r"""NASA polynomials written as :math:`A \times cp + B \times H + C \times S` with,

    - :math:`c_p` heat capacity,
    - :math:`H` enthalpy,
    - :math:`S` entropy.

    Parameters
    ----------
    temperature :
        Temperature
    a1 :
        First NASA coefficients for the polynomials :math:`a_1`
    a2 :
        Second NASA coefficients for the polynomials :math:`a_2`
    a3 :
        Third NASA coefficients for the polynomials :math:`a_3`
    a4 :
        Fourth NASA coefficients for the polynomials :math:`a_4`
    a5 :
        Fifth NASA coefficients for the polynomials :math:`a_5`
    a6 :
        Sixth NASA coefficients for the polynomials :math:`a_6`
    a7 :
        Seventh NASA coefficients for the polynomials :math:`a_7`
    quantity :
        thermodynamic quantity to fit


    Returns
    -------
     : float
        quantity to be computed



    """
    R = ct.gas_constant
    if quantity == 'cp':
        A = 1
        B = 0
        C = 0

    elif quantity == 'enthalpy':
        A = 0
        B = 1
        C = 0

    elif quantity == 'entropy':
        A = 0
        B = 0
        C = 1

    else:
        logger.error('ERROR ! No valid keyword for NASA function')
        quit()

    return R * \
           (A * (a1 +
                 a2 * temperature +
                 a3 * temperature ** 2 +
                 a4 * temperature ** 3 +
                 a5 * temperature ** 4)
            + B * (a1 +
                   a2 * temperature / 2 +
                   a3 * temperature ** 2 / 3 +
                   a4 * temperature ** 3 / 4 +
                   a5 * temperature ** 4 / 5 +
                   a6 / temperature) * temperature
            + C * (a1 * np.log(temperature) +
                   a2 * temperature +
                   a3 * temperature ** 2 / 2 +
                   a4 * temperature ** 3 / 3 +
                   a5 * temperature ** 4 / 4 +
                   a7))


def NASA_leastsq_function(params, *args):
    """Least square method to estimate the best parameters from a given function to fit a set of data

    Parameters
    ----------
    params :
        contains the guessed NASA coefficients for the polynomials (estimated from one of the species)
    *args :
        Contains Temperature, Values of cp, H and S in a single array,
        the A, B and C coefficients from the NASA function

    Returns
    -------
     : float
        value to be minimize for the fit of NASA coefficient


    """
    a1 = params[0]
    a2 = params[1]
    a3 = params[2]
    a4 = params[3]
    a5 = params[4]
    a6 = params[5]
    a7 = params[6]
    T = args[0]
    y = args[1]
    n1 = int(len(T) / 3)
    n2 = 2 * n1

    yfit = np.empty(T.shape)
    yfit[:n1] = NASA(T[:n1], a1, a2, a3, a4, a5, a6, a7, 'cp')
    yfit[n1:n2] = NASA(T[n1:n2], a1, a2, a3, a4, a5, a6, a7, 'enthalpy')
    yfit[n2:] = NASA(T[n2:], a1, a2, a3, a4, a5, a6, a7, 'entropy')

    return y - yfit


def fit_NASA(mechanism, composition_dict={}, spec_arrhenius_dict={}):
    """Function that fits the NASA coefficients for a dictionary of species for a range of temperature

    Parameters
    ----------
    mechanism :
        class `Cantera.mechanisms.Mechanism` object with the species to lump
    composition_dict : dict
        dictionary of species involved in the mixture
        on which we have to determine the NASA coefficient (Default value = {})
    spec_arrhenius_dict : dict
        dictionary of species Arrhenius ratios (Default value = {})

    Returns
    -------
    a_high : list
        list of coefficients for high temperature
    a_low : list
        list of coefficients for low temperature


    """

    ctmech = mechanism.ctmech

    if composition_dict:
        present_species = list(composition_dict.keys())
    else:
        present_species = list(spec_arrhenius_dict.keys())

    index_ref_species = [ctmech.species_index(spec) for spec in present_species]

    # Define the Temperature arrays for the two zones of teh NASA

    Tmin = 300
    Tmax = 5000
    T_rupture = 1000

    list_composition_dict = []

    temperature_high = list(np.linspace(T_rupture, Tmax, 201))
    temperature_low = list(np.linspace(Tmin, T_rupture, 201))

    cp_high = []
    h_high = []
    s_high = []
    cp_low = []
    h_low = []
    s_low = []

    for index_T, T in enumerate(temperature_low):

        if composition_dict:
            list_composition_dict = [composition_dict] * len(temperature_low + temperature_high)

        else:
            composition_dict = {}

            for spec in spec_arrhenius_dict:
                composition_dict[spec] = Arrhenius_function(T, *spec_arrhenius_dict[spec])  # A checker
                list_composition_dict.append(composition_dict)

        cp, h, s = get_thermo(mechanism, T, composition_dict=list_composition_dict[index_T], option='all')

        cp_low.append(cp)
        h_low.append(h)
        s_low.append(s)

    y_low = cp_low + h_low + s_low
    T_low = np.array(temperature_low * 3)

    # Create the variables for the least square method
    for index_T, T in enumerate(temperature_high):
        if composition_dict:
            list_composition_dict = [composition_dict] * len(temperature_low + temperature_high)

        else:
            composition_dict = {}

            for spec in spec_arrhenius_dict:
                composition_dict[spec] = Arrhenius_function(T, *spec_arrhenius_dict[spec])  # A checker
                list_composition_dict.append(composition_dict)

        cp, h, s = get_thermo(mechanism, T, composition_dict=list_composition_dict[index_T + len(temperature_low)],
                              option='all')

        cp_high.append(cp)
        h_high.append(h)
        s_high.append(s)

    y_high = cp_high + h_high + s_high

    T_high = np.array(temperature_high * 3)

    a0_high = []
    a0_low = []
    for i in range(7):
        a0_high.append(ctmech.species()[index_ref_species[0]].thermo.coeffs[i + 1])
        a0_low.append(ctmech.species()[index_ref_species[0]].thermo.coeffs[i + 8])

    params0_high = a0_high
    params0_low = a0_low

    args_high = (T_high, y_high)
    args_low = (T_low, y_low)

    # Lest square method to estimate the NASA coefficients

    result_high = optimize.leastsq(NASA_leastsq_function, params0_high, args=args_high)
    result_low = optimize.leastsq(NASA_leastsq_function, params0_low, args=args_low)

    # Writing of the NASA coefficients

    a_high = []
    a_low = []
    for i in range(7):
        a_high.append(result_high[0][i])
        a_low.append(result_low[0][i])

    # Continuity error calculation on cp
    cp_rupture_high = NASA(1000, *a_high, 'cp')
    cp_rupture_low = NASA(1000, *a_low, 'cp')

    error_Cp = abs(cp_rupture_high - cp_rupture_low) / min(cp_rupture_high, cp_rupture_low)

    logger.debug('Error on Cp at 1000 K with NASA fit = ' + str(error_Cp))

    # Continuity error calculation on enthalpy
    enthalpy_rupture_high = NASA(1000, *a_high, 'enthalpy')
    enthalpy_rupture_low = NASA(1000, *a_low, 'enthalpy')

    error_enthalpy = abs(enthalpy_rupture_high - enthalpy_rupture_low) / min(enthalpy_rupture_high,
                                                                             enthalpy_rupture_low)

    logger.debug('Error on Enthalpy at 1000 K with NASA fit = ' + str(error_enthalpy))

    # Continuity error calculation on entropy
    entropy_rupture_high = NASA(1000, *a_high, 'entropy')
    entropy_rupture_low = NASA(1000, *a_low, 'entropy')

    error_entropy = abs(entropy_rupture_high - entropy_rupture_low) / min(entropy_rupture_high, entropy_rupture_low)

    logger.debug('Error on Entropy at 1000 K with NASA fit = ' + str(error_entropy))

    return a_high, a_low


def compute_rates_profiles(case, mechanism=None):
    """Compute the net reaction rates profiles for each species

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object

    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)

    Returns
    -------
    rates_dictionary : dict
        dictionary containing the net rates of reactions along the case grid for each species

    """
    if not mechanism:
        mechanism = case.mechanism

    rates_dictionary = {}
    for index_reac in range(mechanism.nr):
        rates_dictionary[index_reac] = []

    data_dict = case.data_dict(mechanism)

    grid = data_dict['grid']
    temperature = data_dict['T']
    pressure = data_dict['P']
    mass_fractions = data_dict['Y']

    for index, value in enumerate(grid):

        _, rates = compute_rates(mechanism, temperature[index], pressure[index], mass_fractions[index])

        for index_reac in range(mechanism.nr):
            rates_dictionary[index_reac].append(rates[index_reac])

    return rates_dictionary


def compute_rates(mechanism, temperature, pressure, mass_fractions_list):
    """Retrieves the net reaction rates for a specific thermodynamic state
    Only way to retrieve them for a mechanism with QSS species

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object to be used
    temperature :
        temperature of the gas
    pressure :
        pressure of the gas
    mass_fractions_list :
        list of species mass fractions

    Returns
    -------
    net_production_rates : ndarray
        array of the net production rates
    rates_of_progress :
        array of the rates of progress


    """

    if mechanism.nqss > 0:
        if not hasattr(mechanism, 'skeletal') or not mechanism.skeletal:
            logger.error('ERROR! This function requires the mechanism object to have a skeletal mechanism')
            quit()

        skeletal_ctmech = mechanism.skeletal.ctmech
        ctmech = mechanism.ctmech

        # Reconstructing the composition of the skeletal mechanism
        mass_fractions = {}
        for index, spec in enumerate(ctmech.species_names):
            mass_fractions[spec] = mass_fractions_list[index]

        # Setting the skeletal solution object to correct state
        skeletal_ctmech.TPY = temperature, pressure, mass_fractions
        ctmech.TPY = temperature, pressure, mass_fractions_list

        # Retrieving the concentrations
        concentrations = ctmech.concentrations

        skeletal_concentrations = {}
        for index, spec in enumerate(ctmech.species_names):
            skeletal_concentrations[spec] = concentrations[index]

        # Retrieving the rate constants
        f_rate_constants = list(skeletal_ctmech.forward_rate_constants)
        r_rate_constants = list(skeletal_ctmech.reverse_rate_constants)

        rate_constants = f_rate_constants + r_rate_constants
        # for index, rate in enumerate(r_rate_constants):
        #     if skeletal_ctmech.is_reversible(index):
        #         rate_constants.append(rate)

        # Retrieving the third body values
        third_bodies = []
        for index, reac in enumerate(skeletal_ctmech.reactions()):
            if skeletal_ctmech.reaction_type(index) in [2, 4]:
                efficiencies = reac.efficiencies
                if efficiencies:
                    third_body_list = [(efficiencies[spec] - 1) * skeletal_concentrations[spec]
                                       for spec in list(efficiencies.keys()) if efficiencies[spec] != 0]
                    third_body_list.append(sum(list(skeletal_concentrations.values())))
                    third_body_sum = sum(list(skeletal_concentrations.values()))
                else:
                    third_body_sum = sum(list(skeletal_concentrations.values()))

                third_bodies.append(third_body_sum)
            else:
                third_bodies.append(0)

        # Extracting the QSS routine from the f90 file
        qss_concentrations = compute_qss_concentrations(mechanism, skeletal_concentrations,
                                                        third_bodies, rate_constants)

        for index, spec in enumerate(mechanism.species_qss_names):
            skeletal_concentrations[spec] = qss_concentrations[index] / 1000

        concentrations_list = [skeletal_concentrations[spec] for spec in skeletal_ctmech.species_names]

        skeletal_ctmech.TP = temperature, pressure
        skeletal_ctmech.concentrations = concentrations_list
        rates_of_progress = skeletal_ctmech.net_rates_of_progress
        net_production_rates = skeletal_ctmech.net_production_rates

    else:

        if mechanism.f90:
            ctmech = mechanism.skeletal.ctmech
        else:
            ctmech = mechanism.ctmech

        ctmech.TPY = temperature, pressure, mass_fractions_list
        rates_of_progress = ctmech.net_rates_of_progress
        net_production_rates = ctmech.net_production_rates

    return net_production_rates, rates_of_progress


def compute_qss_concentrations(mechanism, concentrations, third_bodies, rate_constants):
    """Computes the concentrations of the QSS species in an ARC mechanism

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    concentrations :
        concentrations of the transported species
    third_bodies :
        values of the third bodies coefficients
    rate_constants :
        rate constants

    Returns
    -------
    qss_concentrations : ndarray
        concentrations of the QSS species


    """

    ctmech = mechanism.skeletal.ctmech

    # Re-create linear system for QSS species
    species_qss_names = mechanism.species_qss_names
    species_names_full = ctmech.species_names

    ns = mechanism.ns
    nr = mechanism.nr
    nqss = mechanism.nqss

    ns_skel = mechanism.skeletal.ns
    nr_skel = mechanism.skeletal.nr

    import ARCANE.networks as networks
    net = networks.Network(mechanism.skeletal)

    nur = net.nur
    nup = net.nup

    # isqss: array of size ns, True if corresponding spec in species is QSS
    isqss = np.isin(species_names_full, species_qss_names)

    # qss_index_global: array of size nqss, containing index of qss species in global species array
    qss_index_global = np.array(np.argwhere(isqss == True))
    qss_index_global = [value[0] for value in qss_index_global]

    # Get connectivity between species and reactions
    deltass_qss = net.deltaSS
    deltasr_qss = net.deltaSR

    # Reversibility
    is_reversible = []
    for i in range(nr_skel):
        if mechanism.skeletal.ctmech.is_reversible(i):
            is_reversible.append(1)
        else:
            is_reversible.append(0)

    # Trim to only keep QSS/QSS connectivity
    j = 0
    for i in range(ns_skel):
        if i not in qss_index_global:
            deltass_qss = np.delete(deltass_qss, j, 0)
            deltass_qss = np.delete(deltass_qss, j, 1)
            deltasr_qss = np.delete(deltasr_qss, j, 0)
        else:
            j += 1

    # non zero entries in QSS/QSS connectivity matrix
    index_coupled_qss_i, index_coupled_qss_j = np.nonzero(deltass_qss)
    index_qss_reactions_i, index_qss_reactions_j = np.nonzero(deltasr_qss)

    # Retrieving the coupling of the QSS species
    coupled_qss_matrix = np.zeros([nqss, nqss])
    for i in range(nqss):
        for j in index_coupled_qss_j[index_coupled_qss_i == i]:
            if j != i:
                for k in index_qss_reactions_j[index_qss_reactions_i == i]:
                    if ((nup[qss_index_global[i], k] >= 1 and nur[qss_index_global[j], k] >= 1) or
                            ((nur[qss_index_global[i], k] >= 1 and nup[qss_index_global[j], k] >= 1) and
                             is_reversible[k] == 1)):
                        coupled_qss_matrix[i, j] = 1

    # Initialization of the numpy arrays storing the QSS linear system AX = B
    coupled_terms_matrix = np.zeros([nqss, nqss])
    constant_terms_vector = np.zeros(nqss)

    # Connectivity matrices
    deltaReactantsR = abs(nur)
    deltaProductsR = abs(nup)

    j = 0
    for i in range(ns_skel):
        if i not in qss_index_global:
            deltaReactantsR = np.delete(deltaReactantsR, j, 0)
            deltaProductsR = np.delete(deltaProductsR, j, 0)
        else:
            j += 1

    # Removing inert species
    for i in range(deltaReactantsR.shape[0]):
        for j in range(deltaReactantsR.shape[1]):
            if deltaReactantsR[i, j] == - deltaProductsR[i, j]:
                deltaReactantsR[i, j] = 0
                deltaProductsR[i, j] = 0

    deltaReactantsR_qss = deltaReactantsR
    deltaProductsR_qss = deltaProductsR

    indrri, indrrj = np.nonzero(deltaReactantsR_qss)
    indrpi, indrpj = np.nonzero(deltaProductsR_qss)

    for index_spec, qss in enumerate(species_qss_names):

        other_qss_species = species_qss_names.copy()
        other_qss_species.remove(qss)

        denominator = 0
        numerator = 0

        # If matching coefficients
        if indrrj[indrri == index_spec].any() or indrpj[indrpi == index_spec].any():

            # Denominator
            # Forward
            for index_reac in indrrj[indrri == index_spec]:
                reactants = ctmech.reaction(index_reac).reactants
                reactants_keys = list(reactants.keys())

                if not set(other_qss_species).intersection(set(reactants_keys)):

                    reac_rate = rate_constants[index_reac]

                    for spec in reactants_keys:
                        if spec != qss:
                            reac_rate *= concentrations[spec] ** reactants[spec]

                    if ctmech.reaction_type(index_reac) == 2:
                        reac_rate *= third_bodies[index_reac]

                    denominator += reac_rate

            # Backward
            for index_reac in indrpj[indrpi == index_spec]:
                if is_reversible[index_reac] == 1:
                    products = ctmech.reaction(index_reac).products
                    products_keys = list(products.keys())

                    if not set(other_qss_species).intersection(set(products_keys)):

                        reac_rate = rate_constants[index_reac + nr_skel]

                        for spec in products_keys:
                            if spec != qss:
                                reac_rate *= concentrations[spec] ** products[spec]

                        if ctmech.reaction_type(index_reac) == 2:
                            reac_rate *= third_bodies[index_reac]

                        denominator += reac_rate

            # Numerator
            # Forward
            for index_reac in indrpj[indrpi == index_spec]:

                reactants = ctmech.reaction(index_reac).reactants
                reactants_keys = list(reactants.keys())

                if not set(other_qss_species).intersection(set(reactants_keys)):

                    reac_rate = rate_constants[index_reac]

                    for spec in reactants_keys:
                        if spec != qss:
                            reac_rate *= concentrations[spec] ** reactants[spec]

                    if ctmech.reaction_type(index_reac) == 2:
                        reac_rate *= third_bodies[index_reac]

                    numerator += reac_rate

            # Backward
            for index_reac in indrrj[indrri == index_spec]:
                if is_reversible[index_reac] == 1:

                    products = ctmech.reaction(index_reac).products
                    products_keys = list(products.keys())

                    if not set(other_qss_species).intersection(set(products_keys)):

                        reac_rate = rate_constants[index_reac + nr_skel]

                        for spec in products_keys:
                            if spec != qss:
                                reac_rate *= concentrations[spec] ** products[spec]

                        if ctmech.reaction_type(index_reac) == 2:
                            reac_rate *= third_bodies[index_reac]

                        numerator += reac_rate

        if denominator != 0:
            constant_terms_vector[index_spec] = numerator * 1000 / denominator
        else:
            constant_terms_vector[index_spec] = 0

        # Coupled QSS species coefficient
        for j in index_coupled_qss_j[index_coupled_qss_i == index_spec]:
            if coupled_qss_matrix[index_spec, j] == 1 or coupled_qss_matrix[j, index_spec] == 1:

                coefficient = 0
                # Forward
                for k in indrpj[indrpi == index_spec]:
                    reactants = ctmech.reaction(k).reactants
                    reactants_keys = list(reactants.keys())

                    products = ctmech.reaction(k).products

                    if species_qss_names[j] in reactants_keys:
                        reac_rate = - rate_constants[k] * products[species_qss_names[index_spec]]

                        for spec in reactants_keys:

                            if spec not in species_qss_names:
                                reac_rate *= concentrations[spec] ** reactants[spec]

                        if ctmech.reaction_type(index_reac) == 2:
                            reac_rate *= third_bodies[index_reac]

                        coefficient += reac_rate

                # Backward
                for k in indrrj[indrri == index_spec]:
                    if is_reversible[k] == 1:
                        reactants = ctmech.reaction(k).reactants

                        products = ctmech.reaction(k).products
                        products_keys = list(products.keys())

                        if species_qss_names[j] in products_keys:
                            reac_rate = - rate_constants[k + nr_skel] * reactants[species_qss_names[index_spec]]

                            for spec in products_keys:
                                if spec not in species_qss_names:
                                    reac_rate *= concentrations[spec] ** products[spec]

                            if ctmech.reaction_type(index_reac) == 2:
                                reac_rate *= third_bodies[index_reac]

                            coefficient += reac_rate

                coupled_terms_matrix[index_spec, j] = coefficient / (denominator + 1e-60)

    for i, speci in enumerate(species_qss_names):
        coupled_terms_matrix[i, i] = 1

    # Solving the linear system
    qss_concentrations = np.linalg.solve(coupled_terms_matrix, constant_terms_vector)

    return qss_concentrations


def same_species(species, ref_species):
    """Compares two species and returns True if they are the same

    Parameters
    ----------
    species :
        class `Cantera.Species` to be compared
    ref_species :
        class `Cantera.Species` to be compared with

    Returns
    -------
    are_the_same : bool
        True if the species are the same, False of they aren't


    """

    are_the_same = True

    # Comparing composition
    if not species.composition == ref_species.composition:
        are_the_same = False
    else:
        logger.debug("Composition of species " + species.name + ": " + str(species.composition))
        logger.debug("Composition of ref species " + ref_species.name + ": " + str(ref_species.composition))

    if are_the_same:
        # Comparing thermodynamic properties
        nasa_max_difference = 0

        if nasa_max_difference < 1e-10:

            NASA_coeffs = species.thermo.coeffs
            NASA_coeffs_ref = ref_species.thermo.coeffs

            quantities = ['cp', 'enthalpy', 'entropy']
            temperatures = [500, 1000, 2000]

            for quantity in quantities:
                for temperature in temperatures:

                    thermo_value = NASA(temperature, *NASA_coeffs[1:8], quantity)
                    thermo_value_ref = NASA(temperature, *NASA_coeffs_ref[1:8], quantity)

                    thermo_diff = abs((thermo_value - thermo_value_ref)) / abs(thermo_value_ref)
                    logger.debug("Difference on " + quantity + " at " + str(temperature) + " K: " + str(thermo_diff))

                    if thermo_diff > 1:
                        are_the_same = False
                        break

    if are_the_same:
        transport_diffs = []

        if species.transport.geometry != ref_species.transport.geometry:
            are_the_same = False

        # Comparing transport properties
        diameter_diff = (abs(species.transport.diameter
                             - ref_species.transport.diameter)) / (abs(ref_species.transport.diameter) + 1e-60)
        transport_diffs.append(diameter_diff)
        logger.debug("Difference on diameter: " + str(diameter_diff))

        well_depth_diff = (abs(species.transport.well_depth
                               - ref_species.transport.well_depth)) / (abs(ref_species.transport.well_depth) + 1e-60)
        transport_diffs.append(well_depth_diff)
        logger.debug("Difference on well depth: " + str(well_depth_diff))

        rot_relax_diff = (abs(species.transport.rotational_relaxation
                              - ref_species.transport.rotational_relaxation)) / \
                         (abs(ref_species.transport.rotational_relaxation) + 1e-60)
        transport_diffs.append(rot_relax_diff)
        logger.debug("Difference on rotational relaxation: " + str(rot_relax_diff))

        if not all(transport_diffs) < 1:
            are_the_same = False

    return are_the_same


def same_reaction(reaction, ref_reaction):
    """Compares two reactions and returns True if they are the same

    Parameters
    ----------
    reaction :
        class `Cantera.Reaction` to be compared
    ref_reaction :
        class `Cantera.Reaction` to be compared with

    Returns
    -------
    are_the_same : bool
        True if the species are the same, False of they aren't

    """

    are_the_same = True

    # Equation type comparision
    if reaction.reaction_type != ref_reaction.reaction_type:
        are_the_same = False

    # Equation comparision
    if are_the_same:
        if reaction.reactants != ref_reaction.reactants or reaction.products != ref_reaction.products:
            are_the_same = False

        # Comparing Arrhenius coefficients
        if are_the_same:
            if reaction.reaction_type in [1, 2]:
                rate = reaction.rate
                ref_rate = ref_reaction.rate
                A = rate.activation_energy
                b = rate.temperature_exponent
                Ea = rate.activation_energy
                ref_A = ref_rate.activation_energy
                ref_b = ref_rate.temperature_exponent
                ref_Ea = ref_rate.activation_energy

                if [A, b, Ea] != [ref_A, ref_b, ref_Ea]:
                    are_the_same = False

            elif reaction.reaction_type == 4:
                high_rate = reaction.high_rate
                high_ref_rate = ref_reaction.high_rate

                high_A = high_rate.activation_energy
                high_b = high_rate.temperature_exponent
                high_Ea = high_rate.activation_energy
                high_ref_A = high_ref_rate.activation_energy
                high_ref_b = high_ref_rate.temperature_exponent
                high_ref_Ea = high_ref_rate.activation_energy

                low_rate = reaction.low_rate
                low_ref_rate = ref_reaction.low_rate

                low_A = low_rate.activation_energy
                low_b = low_rate.temperature_exponent
                low_Ea = low_rate.activation_energy
                low_ref_A = low_ref_rate.activation_energy
                low_ref_b = low_ref_rate.temperature_exponent
                low_ref_Ea = low_ref_rate.activation_energy

                if [high_A, high_b, high_Ea] != [high_ref_A, high_ref_b, high_ref_Ea] \
                        or [low_A, low_b, low_Ea] != [low_ref_A, low_ref_b, low_ref_Ea]:
                    are_the_same = False

            elif reaction.reaction_type == 5:

                if len(reaction.rates) != len(ref_reaction.rates):
                    are_the_same = False
                else:
                    for index, pressure_rate in enumerate(reaction.rates):

                        pressure = pressure_rate[0]
                        rate = pressure_rate[1]

                        ref_pressure = ref_reaction.rates[index][0]
                        ref_rate = ref_reaction.rates[index][1]

                        A = rate.activation_energy
                        b = rate.temperature_exponent
                        Ea = rate.activation_energy
                        ref_A = ref_rate.activation_energy
                        ref_b = ref_rate.temperature_exponent
                        ref_Ea = ref_rate.activation_energy

                        if [pressure, A, b, Ea] != [ref_pressure, ref_A, ref_b, ref_Ea]:
                            are_the_same = False
                            break

            else:
                logger.warning("Can't assess if the two reactions from the new and th reference mechanism (respectively"
                             + reaction.ID + " and " + ref_reaction.ID + " are the same, they won't be as default.")

    return are_the_same


def create_species(name, initial_species, NASA=None):
    """Creating a new species from an existing one

    Parameters
    ----------
    name : str
        name of the new species
    initial_species : Cantera Species object
        base species to copy properties
    NASA : Cantera NASA object
        user specified NASA polynomial ``ctmech.species()[0]`` object (Default value = None)

    Returns
    -------
    new_species : class `Cantera.Species` object
        new species object from an existing one


    """

    new_species = ct.Species(name, initial_species.composition,
                             initial_species.charge,
                             initial_species.size)

    new_species.transport = initial_species.transport

    if NASA is None:
        new_species.thermo = initial_species.thermo
    else:
        new_species.thermo = NASA

    return new_species


def Zbilger(case, mechanism):
    """Converts a grid array into the :math:`Z` mixture fraction space.

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object used
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
       class :func:`~ARCANE.mechanisms.Mechanism` object used


    Returns
    -------
    z_bilger : `np.ndarray`
        array of Bilger mass fraction


    """

    gas = mechanism.ctmech
    f = ct.CounterflowDiffusionFlame(gas)
    case.restore = 'self'
    case.restore_flame(mechanism, f)

    zc = np.full(f.flame.n_points, 0.0)
    zh = np.full(f.flame.n_points, 0.0)
    zo = np.full(f.flame.n_points, 0.0)

    for j in range(f.flame.n_points):
        f.set_gas_state(j)
        for i in range(gas.n_species):
            zc[j] = gas.n_atoms(gas.species_names[i], 'C') * gas.Y[i] / gas.molecular_weights[i] + zc[j]
            zh[j] = gas.n_atoms(gas.species_names[i], 'H') * gas.Y[i] / gas.molecular_weights[i] + zh[j]
            zo[j] = gas.n_atoms(gas.species_names[i], 'O') * gas.Y[i] / gas.molecular_weights[i] + zo[j]

    # Note : only the transported species are in these calculations
    b = 2 * zc + 0.5 * zh - zo
    z_bilger = np.minimum((b - b[-1]) / (b[0] - b[-1]), np.ones(f.flame.n_points))

    return z_bilger


def extract_scalar(profile, method, grid=None):
    """Extracts a scalar value from an array given a specific method

    Parameters
    ----------
    profile :
        array to treat
    method :
        method to apply, ``'min'``, ``'max'``, ``'init'``, ``'end'``, etc.
    grid :
        (Default value = None)

    Returns
    -------
    scalar : float
        scalar corresponding to the given method


    """

    # Minimum of the curve
    if method in kwdict.methods['min']:
        scalar = np.min(profile)

    # Maximum of the curve
    elif method in kwdict.methods['max']:
        scalar = np.max(profile)

    # Beginning of the curve
    elif method in kwdict.methods['init']:
        scalar = profile[0]

    # End of the curve
    elif method in kwdict.methods['end']:
        scalar = profile[-1]

    # Average of the curve
    elif method in kwdict.methods['mean']:
        scalar = np.mean(profile)

    # Integral of the curve
    elif method in kwdict.methods['int']:
        if grid is None:
            logger.error('A grid is required to extract the integral value')
            quit()
        else:
            scalar = np.trapz(profile, grid)

    else:
        logger.error(f'The specified method ({method}) is not valid, refer to kwdict.methods')

    return scalar


def format_time(seconds):
    """Transforms seconds into human readable string

    Parameters
    ----------
    seconds :
        number of seconds (float)

    Returns
    -------
    string : str
        time into a readable string


    """

    days = int(seconds // 86400)

    hours = int((seconds - days * 86400) // 3600)

    minutes = int((seconds - days * 86400 - hours * 3600) // 60)

    seconds = int(seconds - days * 86400 - hours * 3600 - minutes * 60)

    string = ("{} days, ".format(days) if days else "") + \
             ("{} hours, ".format(hours) if hours else "") + \
             ("{} minutes, ".format(minutes) if minutes else "") + \
             ("{} seconds ".format(seconds) if seconds else "")

    return string


def sort_dict(dictionary, display=True, ascending=True, line_separation=None, line_style='-', format='.6E'):
    """
    Sorts a dictionary by its values

    Parameters
    ----------

    dictionary : dict
        dictionary to sort
    display : bool
        if True, disaplys the ranking
    ascending : bool
        if True, sorts in ascending order else in descending order
    line_separation : float
        value at which to draw a line during display
    line_style : str
        character to be repeated to form the line
    format : str
        values format

    Returns
    -------

    sorted_list : list
        list of tuples sorted according to the given order

    """

    # Selecting the reactions
    sorted_list = sorted(zip(dictionary.values(),
                             dictionary.keys()), key=lambda x: x[0])
    if not ascending:
        sorted_list.reverse()
        line_condition = "dict_value[0] < line_separation"
    else:
        line_condition = "dict_value[0] > line_separation"

    if display:
        length_string = len(max([str(s[1]) for s in sorted_list], key=len))
        line_plotted = False
        for dict_value in sorted_list:
            if line_separation and eval(line_condition) and not line_plotted:
                logger.info(line_style * length_string + line_style * 10)
                line_plotted = True
            logger.info(f"{dict_value[1]:{length_string}} => {dict_value[0]:{format}}")

    return sorted_list
