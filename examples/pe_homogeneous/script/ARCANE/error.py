"""Module implementing error assessment for reduction

All the routine should have the same input and output
"""
import sys

import numpy as np
from scipy.spatial.distance import euclidean

import ARCANE.display as display
import ARCANE.kwdict as kwdict
import ARCANE.tools as tools

logger = display.Logger()
logger.set_log('logReduction')
kwdict = kwdict.Kwdict()


def define_error_quantity(error_type, case, mechanism):
    """Extract the value of a quantity for a case and a mechanism specified

    Parameters
    ----------
    error_type : str
        type of error
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object for comparison

    Returns
    -------
    value : list
        error value for the case and the type of error specified


    Created: 06/12/18 [JW]

    Last modified: 03/01/19 [JW]

    """

    # Initialising extension flags

    if callable(error_type):
        value = error_type(case, mechanism)

    else:

        # Quantities for data
        data_dict = case.data_dict(mechanism)

        # Extracts the method from the quantity string
        error_type_buffer, method, type_of_error, location_of_error = extract_quantity_method(error_type)

        if error_type_buffer:
            error_type = error_type_buffer

        # Updating kwdict with the species_names
        kwdict.update('names', mechanism.kwdict_names_extension)
        kwdict.update('units', mechanism.kwdict_units_extension)

        # Extending data if necessary
        if kwdict.get_base_name(error_type) not in data_dict \
                and error_type not in kwdict.names['tig'] + kwdict.names['sl'] + kwdict.names['thickness']:
            data_dict = case.data_dict(mechanism, extend_data=True)

        # Direct outputs of the vectors
        if kwdict.get_base_name(error_type) in data_dict:
            value = data_dict[kwdict.get_base_name(error_type)]

        # Modified/calculated variables
        elif error_type in kwdict.names['tig']:
            if not method:
                method = 'hr'

            value = tools.get_ignition_delay_time(case, mechanism, method=method)

        elif error_type in kwdict.names['sl']:
            if not method:
                method = 'inlet'

            value = tools.get_laminar_flame_velocity(case, mechanism, method=method)

        elif error_type in kwdict.names['thickness']:
            if not method:
                method = 'thermal'

            value = tools.get_thickness(case, mechanism, method=method)

        else:
            logger.error(f'The quantity {error_type} is not available for that mechanism, error will not be computed')
            value = 'not_valid'

    return value


def error_dictionary(cases_list, ref_mechanism, relevant_var, current_mechanism, max_error_tol, error_factor=1):
    """Calculate the errors between the reference and the current mechanism for all the error specified by the user

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    ref_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        reference class :func:`~ARCANE.mechanisms.Mechanism` object for comparison
    relevant_var : str
        relevant variable ('species', 'reaction' or 'qss')
    current_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        current class :func:`~ARCANE.mechanisms.Mechanism` object for comparison
    max_error_tol : list
        list of tolerances on errors
    error_factor : float
        factor to apply on the error to define the limiting one (Default value = 1)

    Returns
    -------
    error_values : dict
        dictionary of errors
    error_super_list : list
        list of all errors ordered as in cases_list


    Created: 06/12/18 [JW]

    Last modified: 03/01/19 [JW]

    """

    logger.info(str(relevant_var) + " = " + str(getattr(current_mechanism, relevant_var)))
    logger.info('--------------------------')

    # Instantiation
    error_values = {}  # dictionary of errors
    loop = 0
    error_super_list = []
    quantity = ''
    method = ''

    for index_case, case in enumerate(cases_list):
        error_list = case_error(case, ref_mechanism, mechanism=current_mechanism, error_factor=error_factor)

        error_super_list += error_list

        error_dict = case.error_dict

        for index_error, err in enumerate(error_dict):

            if not callable(err):
                # Extracting quantity and method from string
                err_split = err.split(' ')
                last_split = err_split[-1]

                quantity = err
                method = 'none'

                # Checking if last word is a method
                for name in kwdict.methods:
                    if last_split in kwdict.methods[name]:
                        method = last_split
                        quantity = ' '.join(err_split[:-1])

            a = getattr(current_mechanism, relevant_var), case.reactor_type, quantity, \
                method, '{0:2,.2f}'.format(max_error_tol[loop] * 100.)
            loop = loop + 1
            if a not in error_values:
                error_values[a] = []
            error_values[a].append(error_list[index_error])

    logger.info('--------------------------\n')

    return error_values, error_super_list


def cases_list_error(cases_list, reference_mechanism, mechanism=None, error_factor=1, display_results=True,
                     display_data=False):
    """
    Determines the global error on a list of cases

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects
    reference_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        reference class :func:`~ARCANE.mechanisms.Mechanism` object for comparison
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        current class :func:`~ARCANE.mechanisms.Mechanism` object for comparison (Default value = None)
    error_factor : float
        factor to apply on the error to define the limiting one (Default value = 1)
    display_results : bool
        whether results are displayed or not (Default value = True)
    display_data : bool
        whether values used for the error are displayed or not (Default value = False)

    Returns
    -------
    sum_square : float
        global error value
    full_error_dict : dict
        dictionary with all the error values for all cases

    """
    sum_square = 0

    full_error_dict = {}
    for case in cases_list:
        error_storing_dict = case_error(case, reference_mechanism, mechanism=mechanism,
                                        error_factor=error_factor, display_results=display_results,
                                        output_as_dict=True, display_data=display_data)
        full_error_dict[case.myid] = error_storing_dict

        # Computing sum of the squared errors weighed by their tolerances
        weighed_error_list = [(error_storing_dict[key] / case.error_dict[key]) ** 2
                              for key in error_storing_dict]
        sum_square += sum(weighed_error_list)

    return sum_square, full_error_dict


def case_error(ref_case, ref_mechanism, case=None, mechanism=None, error_factor=1,
               display_results=True, output_as_dict=False, display_data=False):
    """Computes the error between 2 mechanisms on error_dict entries

    Parameters
    ----------
    ref_case : class :func:`~ARCANE.cases.Case` object
        reference class :func:`~ARCANE.cases.Case` object for comparison
    ref_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        reference class :func:`~ARCANE.mechanisms.Mechanism` object for comparison
    case : class :func:`~ARCANE.cases.Case` object
        current class :func:`~ARCANE.cases.Case` object for comparison (Default value = None)
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        current class :func:`~ARCANE.mechanisms.Mechanism` object for comparison (Default value = None)
    error_factor : float
        factor to apply on the error to define the limiting one (Default value = 1)
    display_results : bool
        if True, displays results in terminal (Default value = True)
    output_as_dict : bool
        if True, stores the errors as a dictionary instead of list  (Default value = False)`
    display_data : bool
        if True, displays the values used to compute the error

    Returns
    -------
    error_list : list
        list of errors on the case

    """

    if not case and not mechanism:
        logger.error('Error ! Both case and mechanism are empty, specify at least one so that a comparision can be done.')
        quit()

    if not case:
        case = ref_case

    if not mechanism:
        mechanism = ref_mechanism

    error_list = []

    if display_results:
        case_ppt = kwdict.reactor_labels[ref_case.reactor_type][0]
        if hasattr(ref_case, 'pressure'):
            case_ppt += ', p = ' + '{:.2e}'.format(ref_case.pressure)
        if hasattr(ref_case, 'temperature'):
            case_ppt += ', T = ' + '{:.0f}'.format(ref_case.temperature)
        if hasattr(ref_case, 'phi'):
            case_ppt += ', phi = ' + '{:.2f}'.format(ref_case.phi)
        elif hasattr(ref_case, 'composition') and len(ref_case.composition) < 5:
                case_ppt += ', composition = ' + str(ref_case.composition)
        if hasattr(ref_case, 'strain_rate'):
            case_ppt += ', a = ' + '{:.2e}'.format(ref_case.strain_rate)
        if hasattr(ref_case, 'scalar_dissipation_rate'):
            case_ppt += ', chist = ' + '{:.2e}'.format(ref_case.scalar_dissipation_rate)

        if hasattr(ref_case, 'fuel'):
            for f in ref_case.fuel:
                case_ppt = case_ppt + ', ' + f + ' = ' + '{:.0f}'.format(int(ref_case.fuel[f] * 100)) + '%'

        if hasattr(ref_case, 'oxidizer'):
            for o in ref_case.oxidizer:
                case_ppt = case_ppt + ', ' + o + ' = ' + '{:.0f}'.format(int(ref_case.oxidizer[o] * 100)) + '%'

        logger.info('=> on ' + case_ppt)

    error_dict = ref_case.error_dict.copy()

    error_storing_dict = {}
    if ref_case.success:

        for index_error, err in enumerate(error_dict):

            error_val, det, red = compute_error(err, ref_case, ref_mechanism, case=case, mechanism=mechanism)

            if output_as_dict:
                error_storing_dict[err] = error_val
            else:
                error_list.append(abs(error_val))

            if display_results:
                quantity, method, type_of_error, location_of_error = extract_quantity_method(err)

                if method == 'none':
                    method = ''

                # Updating kwdict with the species_names
                kwdict.update('names', mechanism.kwdict_names_extension)
                kwdict.update('units', mechanism.kwdict_units_extension)

                error_name = kwdict.get_full_name(quantity)

                if type_of_error == 'absolute':
                    unit = kwdict.units[error_name].replace('$', '')
                    output_format = '.4e'
                else:
                    error_val *= 100
                    error_dict[err] *= 100
                    unit = '%'
                    output_format = '.5f'

                warning_on = False

                if callable(quantity):
                    error_name = quantity.__name__

                if abs(error_val) > error_dict[err] * error_factor:
                    if error_factor != 1:
                        factor_string = f' * {error_factor}'
                    else:
                        factor_string = ''
                    flag_string = f'  > {error_dict[err]:{output_format}} {unit} {factor_string}  <-- Limiting error'
                    warning_on = True
                else:
                    flag_string = ''

                method = ' ' + method
                text = f'    -> error on {error_name}{method} = {error_val:{output_format}} {unit}{flag_string}'

                if warning_on:
                    logger.warning(text)
                else:
                    logger.info(text)

                if display_data:
                    logger.info(f' ref = {det:.4e}')
                    logger.info(f' new = {red:.4e}')
                else:
                    logger.debug(f' ref = {det:.4e}')
                    logger.debug(f' new = {red:.4e}')

    else:
        error_list = [100] * len(error_dict)
        logger.info('As the case failed or was not computed, the error will not be displayed')

    if output_as_dict:
        return error_storing_dict
    else:
        return error_list


def compute_error(error, ref_case, ref_mechanism, case=None, mechanism=None):
    """Computes the error between two data sets

    Parameters
    ----------
    error : str
        error in the format used by error_dict
    ref_case : class :func:`~ARCANE.cases.Case` object
        reference class :func:`~ARCANE.cases.Case` object for comparison
    ref_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        reference class :func:`~ARCANE.mechanisms.Mechanism` object for comparison
    case : class :func:`~ARCANE.cases.Case` object
        current class :func:`~ARCANE.cases.Case` object for comparison (Default value = None)
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        current class :func:`~ARCANE.mechanisms.Mechanism` object for comparison (Default value = None)

    Returns
    -------
    error_val :
        error value
    det :
        reference value
    red :
        current value

    """

    if not case and not mechanism:
        logger.error('Error ! Both case and mech are empty, specify at least one so that a comparision can be done.')
        quit()

    # Extracting quantity and method
    quantity, method, type_of_error, location_of_error = extract_quantity_method(error)

    if location_of_error:
        real_quantity = quantity
        quantity, method, _, _ = extract_quantity_method(location_of_error)

    if not case:
        case = ref_case

    if not mechanism:
        mechanism = ref_mechanism

    # Take the grid if necessary
    if method in kwdict.methods['int'] + kwdict.methods['dist']:
        x_det = define_error_quantity('grid', ref_case, ref_mechanism)
        x_red = define_error_quantity('grid', case, mechanism)

    # Take the value of the quantity for the reference and the current mechanism
    value_det = define_error_quantity(quantity, ref_case, ref_mechanism)
    value_red = define_error_quantity(quantity, case, mechanism)

    if type(value_det) == str or type(value_red) == str:

        error_val = np.nan
        red = 0
        det = 0

    else:

        # Convert them into arrays
        value_det = np.array(value_det)
        value_red = np.array(value_red)

        # Minimum of the curve
        if method in kwdict.methods['min']:
            det = np.min(value_det)
            red = np.min(value_red)

        # Maximum of the curve
        elif method in kwdict.methods['max']:
            det = np.max(value_det)
            red = np.max(value_red)

        # Beginning of the curve
        elif method in kwdict.methods['init']:
            det = value_det[0]
            red = value_red[0]

        # End of the curve
        elif method in kwdict.methods['end']:
            det = value_det[-1]
            red = value_red[-1]

        # Average of the curve
        elif method in kwdict.methods['mean']:
            det = np.mean(value_det)
            red = np.mean(value_red)

        # Integral of the curve
        elif method in kwdict.methods['int']:
            det = np.trapz(value_det, x_det)
            red = np.trapz(value_red, x_red)

        # Average euclidean distance between the curves (reference is the maximum value)
        elif method in kwdict.methods['dist']:

            try:
                from fastdtw import fastdtw
            except ImportError:
                logger.error('This package is not installed yet. '
                             'Please install it by typing : pip3 install fastdtw.')
                raise

            det = 0
            red, path = fastdtw(value_det, value_red, dist=euclidean)
            ref, path_ref = fastdtw(value_det, value_det*0.0, dist=euclidean)
            error_val = red / np.shape(path)[0] / ref * np.shape(path_ref)[0]

        elif method == 'none':
            det = value_det
            red = value_red

        else:
            logger.error('Error !! Watch out, keyword should be : min, max, init, end, mean, int, dist.')
            sys.exit()

        # Extracting the value at the selected location
        if location_of_error:
            index_of_error_det = list(value_det).index(det)
            index_of_error_red = list(value_red).index(red)

            real_value_det = define_error_quantity(real_quantity, ref_case, ref_mechanism)
            real_value_red = define_error_quantity(real_quantity, case, mechanism)

            det = real_value_det[index_of_error_det]
            red = real_value_red[index_of_error_red]

        if red != 0:
            if np.abs(np.log10(red)) < 10 ** (-10):
                logger.warning('WATCH OUT ! One of your chosen error is too small '
                               'and it leads to big variation of the error.')

        # Calculate the values of the errors and store them in a dictionary
        if method not in kwdict.methods['dist']:
            # Relative error
            if type_of_error == 'relative':
                if red != 0 and det != 0:
                    error_val = (red - det) / det
                else:
                    error_val = 1
            # Absolute error
            elif type_of_error == 'absolute':
                error_val = red - det
            else:
                logger.warning(f'WARNING: This type of error is not available (/{type_of_error}).\n'
                               f'You should use a correct one or code it!. '
                               f'For now, falling back to relative (default).')

        if np.isnan(error_val):
            error_val = 1

    return error_val, det, red


def extract_quantity_method(error):
    """Extracts the strings corresponding to the quantity and the method from error string

    Parameters
    ----------
    error : str
        String describing the error

    Returns
    -------
    quantity : str
        quantity string
    method : str
        method string
    type_of_error : str
        error type string
    error_location : str or None (default)
        error location string
    """

    if not callable(error):
        # Identifying the error type to be used
        type_of_error = 'relative'
        if '#' in error:
            # Extracting quantity and method from string
            treated_error = error.replace(' # ', '#').replace(' #', '#').replace('# ', '#')
            error_split = treated_error.split('#')
            error = error_split[0]
            last_split = error_split[-1]

            for name in kwdict.relativity:
                if last_split in kwdict.relativity[name]:
                    type_of_error = name

        # Identifying the error type to be used
        location_of_error = None
        if '@' in error:
            # Extracting quantity and method from string
            treated_error = error.replace(' @ ', '@').replace(' @', '@').replace('@ ', '@')
            error_split = treated_error.split('@')
            error = error_split[0]
            location_of_error = error_split[-1]

        # Extracting quantity and method from string
        err_split = error.split(' ')
        last_split = err_split[-1]

        quantity = error
        method = 'none'

        # Checking if last word is a method
        for name in kwdict.methods:
            if last_split in kwdict.methods[name]:
                method = last_split
                quantity = ' '.join(err_split[:-1])

        for name in kwdict.compute_methods:
            if last_split in kwdict.compute_methods[name]:
                method = last_split
                quantity = ' '.join(err_split[:-1])

    else:
        quantity = error
        method = 'none'
        type_of_error = 'relative'
        location_of_error = None

    return quantity, method, type_of_error, location_of_error

