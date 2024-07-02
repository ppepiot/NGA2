"""Functions for management of chemical samples and databases"""

# Import statements
import numpy as np
import ARCANE.kwdict as kwdict
import ARCANE.display as display

logger = display.Logger()
logger.set_log('logCompute')

kwdict = kwdict.Kwdict()

def create_samples_database(cases_list, mechanism, data_to_sample=['T', 'P', 'Y', 'HR'],
                            integrated=True, n_samples=None, time_normalised=False):
    """Go through all case instances and extract chemical samples from each of them.
    Results are concatenated into common database

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object containing all info about mechanism
    data_to_sample :
        select physical data to sample (Default value = ['T','P','Y','HR'])
    integrated :
        if False, takes instantaneous data (Default value = True)
    n_samples :
        number of samples (Default value = None)
    time_normalised :
        normalised grid interval (Default value = False)
        
        
    Created: 17/11/14 [PP]
        
        
    Last modified: 17/11/14 [PP] 


    """

    data_to_sample_list = data_to_sample.copy()
    # Extending the list when species are taken
    # In mass fractions
    remake_mass_fraction = False
    intersection = list(set(data_to_sample_list) & set(kwdict.names['Y']))
    if intersection:
        data_to_sample_list.remove(intersection[0])
        data_to_sample_list += mechanism.species_names
        remake_mass_fraction = True

    # In mole fractions
    remake_mole_fraction = False
    intersection = list(set(data_to_sample_list) & set(kwdict.names['X']))
    if intersection:
        data_to_sample_list.remove(intersection[0])
        data_to_sample_list += ['X_' + spec for spec in mechanism.species_names]
        remake_mole_fraction = True

    # In concentrations
    remake_concentrations = False
    intersection = list(set(data_to_sample_list) & set(kwdict.names['c']))
    if intersection:
        data_to_sample_list.remove(intersection[0])
        data_to_sample_list += ['c_' + spec for spec in mechanism.species_names]
        remake_concentrations = True

    # Finding if the grid array is present in the list
    intersection = list(set(data_to_sample_list) & set(kwdict.names['grid'])) \
                   + list(set(data_to_sample_list) & set(kwdict.names['Time'])) \
                   + list(set(data_to_sample_list) & set(kwdict.names['Position']))
    # Fixing the name for easier use in the following
    if intersection:
        data_to_sample_list.remove(intersection[0])
    # Putting the grid as the first item
    data_to_sample_list = ['grid', 'grid interval'] + data_to_sample_list

    samples_database = {}
    for data_name in data_to_sample_list:
        samples_database[data_name] = []
    samples_database['case id'] = []

    if time_normalised:
        samples_database['normalised grid interval'] = []

    for case in cases_list:
        case_samples_dict = extract_samples(case, mechanism, data_to_sample_list, integrated, n_samples)

        for data_name in data_to_sample_list:
            samples_database[data_name] += case_samples_dict[data_name]

        if time_normalised:
            max_time = max(case_samples_dict['grid'])
            samples_database['normalised grid interval'] += [value / max_time for value in case_samples_dict['grid interval']]

        # Adding the case id in the data_dict
        samples_database['case id'] += [case.myid] * len(case_samples_dict['grid'])

    n_data = len(samples_database['grid'])
    if remake_mass_fraction or remake_mole_fraction or remake_concentrations:
        all_mass_samples = []
        all_mole_samples = []
        all_concentrations_samples = []
        for index_sample in range(n_data):
            local_mass_samples = []
            local_mole_samples = []
            local_concentrations_samples = []
            for spec in mechanism.species_names:
                if remake_mass_fraction:
                    local_mass_samples.append(samples_database[spec][index_sample])
                if remake_mole_fraction:
                    local_mole_samples.append(samples_database['X_' + spec][index_sample])
                if remake_concentrations:
                    local_concentrations_samples.append(samples_database['c_' + spec][index_sample])

            # Appending the arrays
            all_mass_samples.append(local_mass_samples)
            all_mole_samples.append(local_mole_samples)
            all_concentrations_samples.append(local_concentrations_samples)

        if remake_mass_fraction:
            samples_database['Y'] = all_mass_samples
        if remake_mole_fraction:
            samples_database['X'] = all_mole_samples
        if remake_concentrations:
            samples_database['c'] = all_concentrations_samples

    # Storing the mechanism
    samples_database['mechanism'] = mechanism

    return samples_database


def extract_samples(case, mechanism, data_to_sample_list, integrated=True, n_samples=None):
    """Extracts samples for dumped profiles

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object containing all info about mechanism
    data_to_sample_list :
        List of physical data to sample
    integrated :
        if False, takes instantaneous data (Default value = True)
    n_samples :
        number of samples (Default value = None)
    
        
    Returns
    -------
    collection of samples

    """
    # Initializing samples collection array
    if not hasattr(case, 'ndata'):
        case.run(mechanism)

    data_dict = case.data_dict(mechanism)
    n_data = case.ndata

    # Check in the case of the function being called alone
    if 'grid' and 'grid interval' not in data_to_sample_list:
        data_to_sample_list = ['grid', 'grid interval'] + data_to_sample_list

    # Checking if the data is consistent
    all_species_in_data = True
    for spec in mechanism.species_names:
        if spec not in data_dict:
            all_species_in_data = False

    if not all_species_in_data:
        case.run(mechanism, overwrite=True)
        data_dict = case.data_dict(mechanism)

    # Fill matrix for data processing
    data_dim = len(data_to_sample_list)
    if not n_samples:
        n_samples = case.nsamples
    elif n_samples == 'all' or n_samples > n_data:
        n_samples = n_data - 1

    data = np.zeros((data_dim, n_data))
    for index, data_name in enumerate(data_to_sample_list):
        data[index, :] = data_dict[data_name]

    # Skipping first result line as the initialization values can cause errors
    data_original = data.copy()
    data = data[:, 1:]
    n_data -= 1

    iterator = 0
    sample_number = 0

    if n_data < n_samples:
        n_samples = n_data

    quotient, rest = divmod(n_data, n_samples)

    sampled_data_dict = {}
    for data_name in data_to_sample_list:
        sampled_data_dict[data_name] = []

    # If there is not at least 2 points per sample, the whole vector must be taken
    if integrated and quotient < 2:
        logger.warning(f'The number of samples ({n_samples}) is too high compared to the number of data ({n_data}) '
                       f'for the integration to take place.')
        integrated = False

    if integrated:

        sampling_arrays = np.zeros((data_dim, quotient))

        # Loop through the rest
        for index in range(n_data):

            line = data[:, index]

            # Arrays for integration
            sampling_arrays[:, iterator] = line

            if (iterator + 1) % quotient == 0 and iterator not in [0, 1]:

                grid_array = sampling_arrays[0, :]

                interval = grid_array[-1] - grid_array[0]

                new_grid_point = (grid_array[-1] + grid_array[0])/2
                sampled_data_dict['grid'].append(new_grid_point)
                sampled_data_dict['grid interval'].append(interval)

                for index_name, data_name in enumerate(data_to_sample_list[2:]):
                    new_data_point = np.trapz(sampling_arrays[index_name + 2, :], grid_array) / interval
                    sampled_data_dict[data_name].append(new_data_point)

                # Reset the arrays
                sampling_arrays = np.zeros((data_dim, quotient))
                iterator = -1
                # quit()
                sample_number += 1

            iterator += 1

    else:

        # Loop through the rest
        for index in range(n_data):

            line = data[:, index]
            # Clip very small concentrations
            line = [0.0 if abs(float(item)) < 1e-15 else item for item in line]

            # Arrays for integration
            sampling_arrays = line

            if (index + 1) % quotient == 0:

                new_grid_point = sampling_arrays[0]
                sampled_data_dict['grid'].append(new_grid_point)

                for index_name, data_name in enumerate(data_to_sample_list[2:]):
                    new_data_point = sampling_arrays[index_name + 2]
                    sampled_data_dict[data_name].append(new_data_point)

                sample_number += 1

        interval_array = [data_original[0, 1] - data_original[0, 0]]
        sampled_data_dict['grid interval'] = interval_array + \
                                             [sampled_data_dict['grid'][index + 1] - sampled_data_dict['grid'][index]
                                              for index, grid_point in enumerate(sampled_data_dict['grid'][1:])]

    return sampled_data_dict


def scalars_database(cases_list, mechanism):
    """Go through all case instances and extract chemical scalars from each of them.
    Results are concatenated into common database

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object containing all info about mechanism
        
    
    Returns
    -------
    collection of samples 


    Created: 17/11/14 [PP]
        
        
    Last modified: 17/11/14 [PP]


    """
    import ARCANE.tools as tools

    scalars_db = np.zeros([len(cases_list), mechanism.ctmech.n_species])
    for index_case, case in enumerate(cases_list):
        data_dict = case.data_dict(mechanism)
        for index_spec, spec in enumerate(mechanism.species_names):
            scalars_db[index_case, index_spec] = tools.extract_scalar(data_dict[spec], 'int', grid=data_dict['grid'])

    return scalars_db
