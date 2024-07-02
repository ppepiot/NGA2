import numpy as np
import cantera as ct
import ARCANE.kwdict as kwdict

kwdict = kwdict.Kwdict()

# Functions directly working on the data_dict that post-processes data

def extend_data_dict(data_dict):
    """Adds new fields to the initial data_dict

    Parameters
    ----------
    data_dict : dict
        initial data_dict

    Returns
    -------
    data_dict :  dict
        extended data_dict

    """
    # Loading mechanism
    mechanism = data_dict['mechanism']

    # Pre-processing of the required values
    n_data = len(data_dict['grid'])
    n_reac = mechanism.nr
    n_spec = mechanism.ns
    elements_in_mech = mechanism.elements_names

    # Compiling mechanism if needed
    if mechanism.f90:
        mechanism.compile()

    # Importing Cantera object and network
    ctmech = mechanism.ctmech
    net = mechanism.network

    # Importing species names
    species_names = ctmech.species_names

    new_data_dict = {}
    # Instantiating the arrays for further computations
    # Single value per point
    for element in elements_in_mech:
        if element not in kwdict.elements_names:
            kwdict.elements_names[element] = 'unknown element'
            kwdict.elements_weights = 0
            kwdict.elements_phi_coeff = 0

        new_data_dict[f'{kwdict.elements_names[element]}'] = np.zeros((n_data))
        new_data_dict[f'{kwdict.elements_names[element]} mass'] = np.zeros((n_data))

    new_data_dict['phi'] = np.zeros((n_data))
    new_data_dict['Prandtl'] = np.zeros((n_data))
    new_data_dict['mass'] = np.zeros((n_data))
    new_data_dict['quantity'] = np.zeros((n_data))

    # Value per reaction per point
    new_data_dict['net rates of progress'] = np.zeros((n_reac, n_data))
    new_data_dict['hr reactions'] = np.zeros((n_reac, n_data))

    # Value per species per point
    new_data_dict['mole fractions'] = np.zeros((n_spec, n_data))
    new_data_dict['concentrations'] = np.zeros((n_spec, n_data))
    new_data_dict['net production rates'] = np.zeros((n_spec, n_data))
    new_data_dict['creation rates'] = np.zeros((n_spec, n_data))
    new_data_dict['destruction rates'] = np.zeros((n_spec, n_data))
    new_data_dict['hr species'] = np.zeros((n_spec, n_data))
    new_data_dict['Schmidt'] = np.zeros((n_spec, n_data))
    new_data_dict['Lewis'] = np.zeros((n_spec, n_data))

    # Looping on all data to retrieve properties
    for index_data in range(n_data):
        # Setting the solution state
        ctmech.TPY = data_dict['T'][index_data], data_dict['P'][index_data], data_dict['Y'][index_data]
        new_data_dict['mole fractions'][:, index_data] = ctmech.X
        new_data_dict['concentrations'][:, index_data] = ctmech.concentrations
        new_data_dict['net production rates'][:, index_data] = ctmech.net_production_rates
        new_data_dict['net rates of progress'][:, index_data] = ctmech.net_rates_of_progress
        new_data_dict['creation rates'][:, index_data] = ctmech.creation_rates
        new_data_dict['destruction rates'][:, index_data] = ctmech.destruction_rates
        new_data_dict['hr reactions'][:, index_data] = - ctmech.net_rates_of_progress * ctmech.delta_enthalpy
        if mechanism.transport != 'Transport':
            new_data_dict['Schmidt'][:, index_data] = ctmech.viscosity / (ctmech.mix_diff_coeffs * ctmech.density)
            new_data_dict['Prandtl'][index_data] = ctmech.viscosity * ctmech.cp_mass / ctmech.thermal_conductivity

        quantity_object = ct.Quantity(ctmech, moles=1)  # 1 mole
        new_data_dict['quantity'][index_data] = quantity_object.mass

        quantity_object = ct.Quantity(ctmech, mass=1)
        new_data_dict['quantity'][index_data] = quantity_object.moles

    if mechanism.transport != 'Transport':
        new_data_dict['Lewis'] = new_data_dict['Schmidt'] / new_data_dict['Prandtl']

    weight = ctmech.molecular_weights

    for index_spec, spec in enumerate(species_names):

        #### HR species
        booli = net.indr[net.inds == index_spec] # reactions containing species i
        for index_reac in booli:
            new_data_dict['hr species'][index_spec, :] += new_data_dict['hr reactions'][index_reac, :]

        ##### Elements profiles computation
        spec_composition = {}
        for element in elements_in_mech:
            spec_composition[element] = ctmech.n_atoms(spec, element)
            ### Number of moles of elements
            new_data_dict[f'{kwdict.elements_names[element]}'] += np.array(data_dict[spec]) \
                                                                  * spec_composition[element] / weight[index_spec]

        #### Mole fractions splitting
        new_data_dict[f'X_{spec}'] = new_data_dict['mole fractions'][index_spec, :]

        #### Concentrations splitting
        new_data_dict[f'c_{spec}'] = new_data_dict['concentrations'][index_spec, :]

        ### Net production rates splitting
        new_data_dict[f'cdot_{spec}'] = new_data_dict['net production rates'][index_spec, :]

        ### Schmidt and Lewis splitting
        new_data_dict[f'Sch_{spec}'] = new_data_dict['Schmidt'][index_spec, :]
        new_data_dict[f'Le_{spec}'] = new_data_dict['Lewis'][index_spec, :]

    # Atom flow_rate
    if 'Velocity' in data_dict:
        for element in elements_in_mech:
            new_data_dict[f'{kwdict.elements_names[element]} flow rate'] = \
                new_data_dict[f'{kwdict.elements_names[element]}'] * data_dict['Velocity']

    new_data_dict['mass'] = 0
    phi_numerator = 0
    phi_denominator = 0
    for element in elements_in_mech:
        new_data_dict[f'{kwdict.elements_names[element]} mass'] = new_data_dict[f'{kwdict.elements_names[element]}'] \
                                                                  * kwdict.elements_weights[element]

        ### Total mass computation
        new_data_dict['mass'] += new_data_dict[f'{kwdict.elements_names[element]} mass']

        # Numerator of phi computation
        if kwdict.elements_phi_coeff[element] > 0:
            phi_numerator += kwdict.elements_phi_coeff[element] * new_data_dict[f'{kwdict.elements_names[element]}']
        else:
            phi_denominator -= kwdict.elements_phi_coeff[element] * new_data_dict[f'{kwdict.elements_names[element]}']

    ### Equivalence ratio computation
    # This first formula uses an atomic conservation formula of the equivalence ratio equivalent to the classical one
    # for CxHy hydrocarbons but yielding different results for oxygenated fuels CxHyOz.
    # new_data_dict['phi'] = (0.5 * new_data_dict['hydrogen'] + 2 * new_data_dict['carbon']) / new_data_dict['oxygen']
    new_data_dict['phi'] = phi_numerator / phi_denominator

    # This one yields the same results as the classical definition by assuming oxidizer as air
    # and oxygen and nitrogen diffusion being the same
    # The sources must be modified to be used
    # new_data_dict['phi'] = (0.5 * new_data_dict['hydrogen'] + 2 * new_data_dict['carbon']
    #                         - new_data_dict['oxygen'] + new_data_dict['nitrogen'] / 3.76) \
    #                        / (new_data_dict['nitrogen'] / 3.76)

    # Updating initial dictionary
    data_dict.update(new_data_dict)

    return data_dict
