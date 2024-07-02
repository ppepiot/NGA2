"""
Module for handling the files and the database

Author: QC 10/06/2022
"""

# Import statements
import subprocess
from pathlib import Path
import os
import shutil
import sys
import numpy as np
import csv
import h5py

import ARCANE.cases as cases
import ARCANE.mechanisms as mechanisms
import ARCANE.display as display
import ARCANE.kwdict as kwdict

logger = display.Logger()
logger.set_log('logCompute')

database_path = '.'
database_system = 'database'
storage_format = 'hdf5' #'ascii'

kwdict = kwdict.Kwdict()


def create_dir(myfolder, overwrite=False):
    """
    Clean up and create folder where all simulations will be stored

    # Check if directory "folder" exists, if not create it.
    # If it exists, stop, unless optional flag allows overwriting
    # Note: if "folder" is a file rather than a folder, save it with a different name first.

    Parameters
    ----------
    myfolder : str
        Name of folder to store all simulations
    overwrite : bool
        if True, existing folder will be overwritten

    Returns
    -------
    success : True or None
        True if successful creation


    Created: 11/12/17 [PP]
    Last modified: 23/06/23 [PP]

    """

    # TODO - Write unit test for function create_dir
    # TODO - can we use only os instead of os and pathlib? (or only pathlib)

    if Path(myfolder).exists() and not Path(myfolder).is_dir():
        logger.info("Warning creating case database: ', myfolder, ' is a file, not a folder! ")
        logger.info("... Moved to ', myfolder.strip(), '_old")
        movefile = myfolder + '_old'
        os.rename(myfolder, movefile)

    try:
        if not Path(myfolder).exists():
            # Folder does not exist, just create it
            os.mkdir(myfolder)
            return True
        elif Path(myfolder).exists() and overwrite is True:
            # Folder exists, but should be overwritten
            shutil.rmtree(myfolder)
            os.mkdir(myfolder)
            return True
        else:
            # If folder already exists and is not to be overwritten
            raise TypeError

    except TypeError:
        logger.debug("Folder {0} already exists and will not be overwritten: " \
                     "data will be re-used.".format(myfolder))


def convert_ascii_to_hdf5(mechanism, directory='./new_database'):
    """
    Converts all the solutions of a mechanism into hdf5 storage format
    By default, it will create in the current directory a 'new_database' directory
    in which the new mechanisms directories are stored.
    You will have to manually change its name to 'database' in order to use it.

    Parameters
    ----------
    mechanism: class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object that will be converted
    directory : str
        name of the new database directory

    """
    global storage_format
    if storage_format != 'ascii':
        storage_format = 'ascii'

    # Retrieving the list of the profile files
    cases_list_profiles = mechanism.computed_cases(display=False, solution=False, create_case=True)

    # Creation of the new database folder
    create_dir(directory)

    global database_path
    old_database = database_path
    database_path = directory

    storage_format = 'hdf5'

    old_mechanism_directory = mechanism.mechdir

    # Create new mechanism directory
    create_dir(f"{database_path}/{mechanism.name}")

    # Writing the mechanism files
    mechanism.write(f"{database_path}/{mechanism.name}/{mechanism.name}.cti")
    if mechanism.f90:
        mechanism.write_f90(f"{database_path}/{mechanism.name}/{mechanism.name}.f90")

    # Converting the profile files
    for case in cases_list_profiles:
        # Load and read ascii file
        mechanism.mechdir = old_mechanism_directory
        data_dict = read_ascii_profiles(case, mechanism)
        # Write the hdf5 file
        mechanism.mechdir = f"{database_path}/{mechanism.name}"
        write_hdf5_profiles(data_dict, case, mechanism)

        if case.myid.startswith('C1D'):
            convert_xml_solution(case, mechanism=mechanism, old_database=old_database, new_database=directory)


def convert_xml_solution(case, mechanism=None, old_database='database', new_database='new_database'):
    """
    Restores a Flame object from the profile file.

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
    old_database: str
        name of the database from which the solution will be extracted
    new_database: str
        name of the database towards which the solution will be written

    """
    import cantera as ct

    if not mechanism:
        mechanism = case.mechanism

    ctmech = mechanism.ctmech

    if case.reactor_type == 'C1DP':
        flame_object = ct.FreeFlame(ctmech)

    elif case.reactor_type == 'C1DCF':
        flame_object = ct.CounterflowDiffusionFlame(ctmech)

    elif case.reactor_type == 'C1DF':
        flame_object = ct.Flamelet(ctmech)

    elif case.reactor_type == 'C1DB':
        flame_object = ct.BurnerFlame(ctmech)

    else:
        logger.error('This reactor does not exist. You can only rewrite solutions for \
                premixed flames, counterflow diffusion flames, flamelets \
                and burners.')

    loglevel = 0

    load_solutions_names(case, mechanism=mechanism, file_format='ascii')
    solution_path, solution_name = xml_solution_file_and_name(case.xml_solution, file_format='ascii')
    solution_path = solution_path.replace(new_database, old_database)
    flame_object.restore(solution_path, name=solution_name, loglevel=loglevel)

    load_solutions_names(case, mechanism=mechanism, file_format='hdf5')
    solution_path, solution_name = xml_solution_file_and_name(case.xml_solution, file_format='hdf5')
    flame_object.save(solution_path, name=solution_name, loglevel=loglevel)


def get_cases_ids(mechanism, solution=False):
    """
    Retrieves the ids of the computed case for the input mechanism

    Parameters
    ----------
    mechanism: class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object that will be converted
    solution: bool
        if True, check the files sol_***.xml and not the files profile_***

    Returns
    -------
    computed_cases_ids : list
        list of ids of the computed cases

    """

    if database_system == 'database':
        if storage_format == 'ascii':
            dir_content = os.walk(mechanism.mechdir)
            dir_content = [val for val in dir_content]
            files = dir_content[0][2]
            if solution:
                computed_cases_id = [file for file in files if file.startswith('sol_')]
                computed_cases_id = [case_str.replace('sol_', '') for case_str in computed_cases_id]
                computed_cases_id = [case_str.replace('.xml', '') for case_str in computed_cases_id]
            else:
                computed_cases_id = [file for file in files if file.startswith('profile_')]
                computed_cases_id = [case_str.replace('profile_', '') for case_str in computed_cases_id]
        else:
            solution_file = f"{mechanism.mechdir}/{mechanism.name}.h5"
            if os.path.isfile(solution_file):
                h5_file = h5py.File(solution_file, 'r+')
                computed_cases_id = list(h5_file.keys())
            else:
                computed_cases_id = []
    else:
        computed_cases_id = os.walk(cases.Case.casedir)
        computed_cases_id = [computed_case_id[0]
                             for computed_case_id in computed_cases_id]
        computed_cases_id = computed_cases_id[1:]

    return computed_cases_id


def load_solutions_names(case, mechanism=None, system=None, file_format=None):
    """
    Instantiating solutions names

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object
    system : str ('database' or 'cases')
        type of data system
    - cases and mechanisms stored in different directories (cases)
    - cases stored inside mechanisms directories (database)
    file_format : str ('ascii' or 'hdf5')
        file format used for storing data

    """

    if not mechanism:
        mechanism = case.mechanism

    if not system:
        system = database_system

    if not file_format:
        file_format = storage_format

    # Instantiating solution name
    if system == 'database':
        if file_format == 'ascii':
            case.solution = f"{mechanism.mechdir}/profile_{case.myid}"
            case.xml_solution = f"{mechanism.mechdir}/sol_{case.myid}.xml"
        else:
            case.solution = f"{mechanism.mechdir}/{mechanism.name}.h5/{case.myid}"
            case.xml_solution = f"{mechanism.mechdir}/{mechanism.name}.xml/{case.myid}"
    else:
        create_dir(f"{case.casedir}/{case.myid}")
        case.solution = f"{case.casedir}/{case.myid}/profile_{mechanism.name}"
        case.xml_solution = f"{case.casedir}/{case.myid}/sol_{mechanism.name}.xml"


def xml_solution_file_and_name(xml_solution_path, file_format=None):
    """
    Parse a xml solution path formatted either for ascci or hdf5 storing
    and returns the path to the file and the name

    Parameters
    ----------
    xml_solution_path : str
        Absolute or relative path of the xml solution for 1D flames with intern solution name
    file_format : str ('ascii' or 'hdf5')
        file format used for storing data


    Returns
    -------
    solution_path : str
        Path to the actual file
    solution_name : str or None
        Name of the intern solution ('hdf5' storing system)

    """

    if not file_format:
        file_format = storage_format

    if file_format == 'hdf5':
        # Separating the group name from the solution file
        splitted_name = xml_solution_path.split('/')
        solution_name = splitted_name[-1]
        solution_path = xml_solution_path.replace(f'/{solution_name}', '')
    else:
        solution_path = xml_solution_path
        solution_name = 'solution'

    return solution_path, solution_name


def check_solution_existence(case, mechanism=None, file_format=None):
    """
    Checks if the solution file exists

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object
    file_format : str ('ascii' or 'hdf5')
        file format used for storing data

    Returns
    -------
    exists: bool

    """

    if not mechanism:
        mechanism = case.mechanism

    if not file_format:
        file_format = storage_format

    load_solutions_names(case, mechanism=mechanism)

    if file_format == 'ascii':
        exists = os.path.isfile(case.solution)
        if case.reactor_type.startswith('C1D'):
            exists = check_xml_solution_existence(case, mechanism=mechanism)
    else:
        # Separating the group name from the solution file
        splitted_name = case.solution.split('/')
        group_name = splitted_name[-1]
        solution_name = case.solution.replace(f'/{group_name}', '')

        if os.path.isfile(solution_name) and group_name in h5py.File(solution_name, 'r+'):
            exists = True
        else:
            exists = False

    return exists

def check_xml_solution_existence(case, mechanism=None):
    """
    Checks if the xml solution file exists

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    exists: bool

    """

    if not mechanism:
        mechanism = case.mechanism

    load_solutions_names(case, mechanism=mechanism)

    exists = os.path.isfile(case.xml_solution)

    return exists

def write_profiles(data_dict, case, mechanism):
    """
    Calls the write function associated with the wanted data type

    Parameters
    ----------
    data_dict : dict
        dictionary of simulation data
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    """

    getattr(sys.modules[__name__], f'write_{storage_format}_profiles')(data_dict, case, mechanism)


def write_ascii_profiles(data_dict, case, mechanism):
    """
    Writes data as an ascii file

    Parameters
    ----------
    data_dict : dict
        dictionary of simulation data
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    """

    # Instantiating solution name
    load_solutions_names(case, mechanism, file_format='ascii')

    if sys.version_info[0] < 3:
        infile = open(case.solution, 'w')
    else:
        infile = open(case.solution, 'w', newline='', encoding='utf8')

    with infile as outfile:
        writer = csv.writer(outfile, delimiter=' ')
        writer.writerow(list(data_dict.keys()))
        for index in range(case.ndata):
            writer.writerow([data_dict[key][index] for key in list(data_dict.keys())])
        if case.ndata == 0:
            writer.writerow(np.zeros(len(list(data_dict.keys()))))
            writer.writerow(np.zeros(len(list(data_dict.keys()))))


def write_hdf5_profiles(data_dict, case, mechanism):
    """
    Writes data as an hdf5 file

    Parameters
    ----------
    data_dict : dict
        dictionary of simulation data
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    """

    # Instantiating solution name
    load_solutions_names(case, mechanism, file_format='hdf5')

    # Separating the group name from the solution file
    splitted_name = case.solution.split('/')
    group_name = splitted_name[-1]
    solution_name = case.solution.replace(f'/{group_name}', '')

    # Opening the file
    if os.path.isfile(solution_name):
        h5_solution = h5py.File(solution_name, 'r+')
    else:
        h5_solution = h5py.File(solution_name, 'w')

    # Creating a group with the case name
    if group_name in h5_solution:
        del h5_solution[group_name]

    h5_solution.create_group(group_name)

    add = h5_solution[group_name]

    # Adding the data to the group
    for variable_name in data_dict:
        add.create_dataset(variable_name, data=data_dict[variable_name])

    h5_solution.close()


def read_profiles(case, mechanism):
    """
    Calls the read function associated with the wanted data type

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    data_dict : dict
        dictionary of simulation data

    """

    data_dict = getattr(sys.modules[__name__], f'read_{storage_format}_profiles')(case, mechanism)

    return data_dict


def read_ascii_profiles(case, mechanism):
    """
    Reads data stored as an ascii file
    
    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    data_dict : dict
        dictionary of simulation data
    """

    # Instantiating solution name
    load_solutions_names(case, mechanism, file_format='ascii')

    # Go through profile file and extract nsamples states (class attribute)
    if sys.version_info[0] < 3:
        infile = open(case.solution, 'rb')
    else:
        infile = open(case.solution, 'r', newline='', encoding='utf8')

    with infile as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')
        line = next(reader)
        if line[0] == '#':
            line = next(reader)

        data_names = line

        data = np.loadtxt(csvfile, delimiter=' ')  # , quotechar='|')

    data_dict = {}
    for index, name in enumerate(data_names):
        data_dict[name] = data[:, index]
        data_dict[kwdict.get_base_name(name)] = data[:, index]

    data_dict['Grid'] = data[:, 0]
    data_dict['grid'] = data[:, 0]

    case.ndata = len(data_dict["grid"])

    return data_dict


def read_hdf5_profiles(case, mechanism):
    """
    Reads data stored in a hdf5 file

    Parameters
    ----------
    case : class :func:`~ARCANE.cases.Case` object
        class :func:`~ARCANE.cases.Case` object
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    data_dict : dict
        dictionary of simulation data
    """

    # Instantiating solution name
    load_solutions_names(case, mechanism, file_format='hdf5')

    # Separating the group name from the solution file
    splitted_name = case.solution.split('/')
    group_name = splitted_name[-1]
    solution_name = case.solution.replace(f'/{group_name}', '')

    # Opening the file
    h5_solution = h5py.File(solution_name, 'r+')

    data_dict = {}
    dict_names = list(h5_solution[group_name].keys())

    for name in dict_names:
        data_dict[name] = h5_solution[group_name][name][()]
        data_dict[kwdict.get_base_name(name)] = h5_solution[group_name][name][()]

    data_dict['Grid'] = data_dict[kwdict.reactor_grid_name[case.reactor_type]]
    data_dict['grid'] = data_dict['Grid']

    case.ndata = len(data_dict["grid"])

    return data_dict


def init_database(path_to_database='.', system='database', file_format='ascii'):
    """
    Set the directory in which all cases and mechanisms will be stored
    This is useful if you want to make a database with all your cases avoiding repetition

    Parameters
    ----------
    path_to_database : str
        path to the database directory
    (you can make it start with $ to retrieve an environment variable)
    system : str ('database' or 'cases')
        type of data system
    - cases and mechanisms stored in different directories (cases)
    - cases stored inside mechanisms directories (database)
    file_format : str ('ascii' or 'hdf5')
        file format used for storing data

    """

    if path_to_database.startswith('$'):
        if path_to_database[1:] in os.environ:
            path_to_database = os.environ.get(path_to_database[1:])
        if not path_to_database:
            logger.error('The environment variable does not exists. Local directory will be used.')
            path_to_database = '.'

    if not os.path.isdir(path_to_database):
        create_dir(path_to_database)

    if system == 'database':
        logger.info('Database storing system will be used.')
        mechanisms.init_mechanism_database(path_to_database + '/database')
        path_to_database += '/database'
    else:
        logger.info('Classical storing system will be used.')
        cases.init_case_database(path_to_database + '/cases')
        mechanisms.init_mechanism_database(path_to_database + '/mechs')

    logger.info(f'Storing directory has been set. All cases and mechanisms data will be stored here: {path_to_database}'
                f' as {file_format} files.\n')

    global database_path
    global database_system
    global storage_format
    database_path = path_to_database
    database_system = system
    storage_format = file_format

