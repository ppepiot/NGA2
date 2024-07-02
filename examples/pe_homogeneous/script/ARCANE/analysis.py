"""Functions for analysis of the chemistry"""

# Import statements

# General imports
import os
import sys

# ARCANE imports
import ARCANE.custom.custom_kinetics as custom_kinetics
import ARCANE.display as display
import ARCANE.drgep as drg
import ARCANE.error as error
import ARCANE.database as database
import ARCANE.lumping as lumping
import ARCANE.mechanisms as mechanisms
import ARCANE.networks as nk
import ARCANE.qss_species as qss
import ARCANE.sampling as sampling
import ARCANE.tools as tools

# Useful to sort a dictionary
import collections
import operator

# Useful to repeat list
from itertools import repeat

# Math imports
import numpy as np
import cantera as ct

# Graph imports
import pandas as pd
import matplotlib.pyplot as plt
from graphviz import Digraph

# Image storing
from PIL import Image

logger = display.Logger()
logger.set_log('logAnalysis')
logger.create_log('logAnalysis')

global analysis_dir
analysis_dir = 'analysis_results'
database.create_dir(analysis_dir, overwrite=False)


def mechanism_analysis(mechanism_list, fuel=None, oxidizer='air', temperature=300,
                       pressure=101325.0, phi=1.0, reactor='HP', plot=True):
    """ Function give features about the mechanism used and a graph with the
    species in the post flame zone.

    :param mechanism_list: list of mechanism objects
    :param fuel: list of fuels used to compute
    :param oxidizer: oxidizer used, default is air
    :param temperature: Temperature of cold gases
    :param pressure: Pressure of the reactor
    :param phi: Equivalence ratio
    :param reactor: Type of reactor. 'HP' for constant pressure, 'UV' for constant volume
    :param plot: Logical to plot hot gases composition histogram

    :return None

    Created: 29/08/2019 [JW]
    Last modification 03/09/2019 [AP]

    """

    if fuel is None:
        logger.error('Error : The fuel list is empty or incomplete.')

    logger.info('---------------------------')
    logger.info('Starting mechanism analysis')
    logger.info('---------------------------')

    for id_mec, mec in enumerate(mechanism_list):
        cold_species = {}
        hot_species = {}
        bars = {}

        logger.info('\n' + str(mec.name))
        logger.info('..................')
        if mec.f90:
            logger.info('This mechanism contains ' + str(mec.ns) + ' species, ' + str(mec.nqss)
                        + ' qss species and ' + str(mec.nr) + ' reactions.')
        else:
            logger.info('This mechanism contains ' + str(mec.ns) + ' species and ' + str(mec.nr) + ' reactions.')

        gas = mec.ctmech
        gas.TP = temperature, pressure
        if oxidizer == 'air':
            gas.set_equivalence_ratio(phi, fuel[id_mec], 'O2:1.0, N2:3.76')
        elif oxidizer == 'oxygen':
            gas.set_equivalence_ratio(phi, fuel[id_mec], 'O2:1.0')
        else:
            logger.error('Error : please, choose air or oxygen.')

        # Information for the cold mixture
        logger.info('For a stoichiometric mixture at ' + str(gas.P) + ' Pa and ' + str(gas.T) + ' K :')
        logger.info('=> The cold mix is :')
        for id_spec, spec in enumerate(gas.species_names):
            cold_species[spec] = gas.X[id_spec]
        cold_species_sorted = sorted(cold_species.items(), key=operator.itemgetter(1), reverse=True)
        cold = collections.OrderedDict(cold_species_sorted)
        for (spec, concentrations) in cold.items():
            if concentrations != 0:
                logger.info('---> ' + str(spec) + ' : ' + '{:.2e}'.format(concentrations))

        # Information for the hot mixture
        gas.equilibrate(reactor)
        for id_spec, spec in enumerate(gas.species_names):
            hot_species[spec] = gas.X[id_spec]
        hot_species_sorted = sorted(hot_species.items(), key=operator.itemgetter(1), reverse=True)
        hot = collections.OrderedDict(hot_species_sorted)

        logger.info('=> The post-combustion temperature is ' + str(gas.T) + ' K')
        logger.info('=> The important species in the post-flame zone are :')
        for (spec, concentrations) in hot.items():
            if concentrations >= 1e-4:
                logger.info('---> ' + str(spec) + ' : ' + '{:.2e}'.format(concentrations))
        bars[mec.name] = list(hot.values())

        if plot:
            # Graph for important species (in log)
            df = pd.DataFrame(bars[mec.name], index=list(hot.keys()))
            df.plot.bar(rot=90, logy=True)
            plt.show()


def difference_mechanism_analysis(mechanism_list, directory=None, nb_line='auto'):
    """ Given two mechanisms, this function creates an HTML file that lists all
    the species that are the same and those that are different, and in which
    mechanism they belong.

    :param mechanism_list: list of mechanism objects
    :param directory: string to make the directory for storing html files
    :param nb_line: choose the number of lines in the table

    :return None

    Created: 2019/08/29 [JW]
    Last modification 2020/01/09 [AP]
      - BUGFIX for empty lst variables

    """

    logger.info('-----------------------------------------------')
    logger.info('Starting difference between mechanisms analysis')
    logger.info('-----------------------------------------------')

    if len(mechanism_list) != 2:
        logger.error("You should have two mechanisms to compare (not less, not more).")
        quit()

    # Create the folder where the html files are gonna be put in
    if not directory:
        directory = analysis_dir + '/analysis_difference_mechanism'
    database.create_dir(directory, overwrite=False)

    analysed_quantities = ['species', 'reactions', 'qss']

    for quantity in analysed_quantities:
        # Get names
        if quantity == 'species':
            list_1 = mechanism_list[0].ctmech.species_names
            list_2 = mechanism_list[1].ctmech.species_names
        elif quantity == 'reactions':
            if mechanism_list[0].nqss > 0:
                list_1 = mechanism_list[0].skeletal.ctmech.reaction_equations()
            else:
                list_1 = mechanism_list[0].ctmech.reaction_equations()
            if mechanism_list[1].nqss > 0:
                list_2 = mechanism_list[1].skeletal.ctmech.reaction_equations()
            else:
                list_2 = mechanism_list[1].ctmech.reaction_equations()
        elif quantity == 'qss':
            if not mechanism_list[0].species_qss_names and not mechanism_list[1].species_qss_names:
                break
            else:
                list_1 = []
                list_2 = []
                if mechanism_list[0].species_qss_names:
                    list_1 = list(mechanism_list[0].species_qss_names)
                if mechanism_list[1].species_qss_names:
                    list_2 = list(mechanism_list[1].species_qss_names)
        else:
            raise ValueError

        logger.info('-> ' + str(quantity) + ' are analysed')

        # Lists that are going to be filled
        common_list = []
        different_list_1 = []
        different_list_2 = []

        # Filling the list with common, different species
        for element_1 in list_1:
            if element_1 in list_2:
                common_list.append(element_1)
        for element_common in common_list:
            list_1.remove(element_common)
            list_2.remove(element_common)
        different_list_1.extend(list_1)
        different_list_2.extend(list_2)

        lists = [common_list, different_list_1, different_list_2]
        list_labels = ['Common ' + str(quantity), 'Different ' + str(quantity) + ' from mechanism 1',
                       'Different ' + str(quantity) + ' from mechanism 2']

        # Number of maximum line in the table
        if nb_line == 'auto':
            nb_line_max = max(len(common_list), len(different_list_1), len(different_list_2))
        else:
            nb_line_max = int(nb_line)

        # Store the data for plotting (table, color and label)
        label_data = []
        for id_label, label in enumerate(list_labels):
            label_data.extend(repeat(label, 1 + (len(lists[id_label])) // (nb_line_max + 1)))

        table_data = [[] * n for n in range(nb_line_max)]
        index = 0
        for id_list, lst in enumerate(lists):
            if not lst:
                lst.append(' ')
            for id_spec, spec in enumerate(lst):
                index = id_spec - nb_line_max * (id_spec // nb_line_max)
                table_data[index].append(spec)
            for id_spec in range(index + 1, nb_line_max, 1):
                table_data[id_spec].append(' ')

        # Print the html file
        show_table(table_data, label_data, directory + '/difference_' + str(quantity) + '.html')


def show_table(table_data, label_data, filename):
    """ Given the table data, the different labels and a filename, build the
    html table for output.

    :param table_data: list of species or reactions
    :param label_data: list of mechanisms to fit the different columns
    :param filename: string

    :return None

    Created: 29/08/2019 [JW]

    """

    table = pd.DataFrame(table_data, columns=label_data)

    path_to_init = tools.__file__
    terminator = path_to_init.rfind('/tools.py')
    path = (path_to_init[:terminator])
    path += '/style_sheets'

    html_string = """
<html>
  <head><title>HTML Pandas Dataframe with CSS</title></head>
  <link rel="stylesheet" type="text/css" href=""" + path + "/table_style.css" + """/>
  <body>
    {table}
  </body>
</html>.
"""

    # Output an html file
    with open(filename, 'w') as f:
        f.write(html_string.format(table=table.to_html(classes='mystyle')))

    logger.info('-> HTML file has been generated under : ' + filename)


def drgep_analysis(case_list, mechanism_list, reduction_type, directory=None, each_case=False,
                   integrity=True, qss_method='LOI', nb_line='auto'):
    """ Given two mechanisms, this function creates an HTML file that shows the
    results of the DRGEP analysis. Each case can be given separately in output
    or together.

    :param case_list: list of case objects
    :param mechanism_list: list of mechanism objects
    :param reduction_type: string type of reduction
    :param directory: string to make the directory for storing html files
    :param each_case: boolean to plot each case or the ensemble of cases
    :param integrity: boolean to choose to analyse integrity or not
    :param nb_line: string or integer corresponding to the number of line in the table.
    Set 'auto' to fit to the number of reaction or species in the mechanism
    :param qss_method: string corresponding to the QSS method: 'LOI' or 'AOI'

    :return None

    Created: 29/08/2019 [JW]
    Last modification: 12/12/2019 [AP]
        - New parameter for the number of line in the table
        - Remove useless lines in case of 'lumping' reduction_type
        - New parameter for the QSS method

    """

    # Create the folder where the html files are gonna be put in
    # Create the folder where the html files are gonna be put in
    if not directory:
        directory = analysis_dir + '/analysis_difference_drgep'
    database.create_dir(directory, overwrite=False)

    # Warning for type of reduction
    if reduction_type == 'species' or reduction_type == 'qss' or reduction_type == 'lumping':
        number_variable = 'ns'
    elif reduction_type == 'reactions':
        number_variable = 'nr'
    else:
        number_variable = ''
        logger.error('Keyword for reduction type should be species or reactions !')

    logger.info("\n=> Building DRGEP for " + reduction_type)

    # Number of maximum line in the table
    if nb_line == 'auto':
        nb_line_max = 0
        for mec in mechanism_list:
            if number_variable == 'ns':
                if nb_line_max < mec.ns:
                    nb_line_max = mec.ns
            elif number_variable < 'nr':
                nb_line_max = mec.nr
    elif type(nb_line) is int:
        nb_line_max = int(nb_line)
    else:
        raise ValueError

    # Data tables that are going to be filled (table, label and color)
    label_data = []
    table_data = [[] * n for n in range(nb_line_max)]
    color_data = [[] * n for n in range(nb_line_max)]

    # Filling the label_data file to reproduce header for different column
    for mec in mechanism_list:
        if each_case:
            for case in case_list:
                if reduction_type == 'lumping':
                    label_data.extend(repeat(mec.name + case.myid,
                                             1 + getattr(mec, number_variable) // (nb_line_max + 1)))
                else:
                    label_data.extend(repeat(mec.name + case.myid,
                                             1 + getattr(mec, number_variable) // (nb_line_max + 1)))
        else:
            if reduction_type == 'lumping':
                label_data.extend(repeat(mec.name, 1 + getattr(mec, number_variable) // (nb_line_max + 1)))
            else:
                label_data.extend(repeat(mec.name, 1 + getattr(mec, number_variable) // (nb_line_max + 1)))

    # Build the different table if drgep is wanted for each case or just
    # overall
    nb_lump = 0
    for id_mec, mec in enumerate(mechanism_list):
        if each_case:
            for case in case_list:
                # Computing is different if the reduction type is different
                # DRGEP for reaction and species, finding isomers for lumping
                # and LOI for qss.
                if reduction_type == 'species' or reduction_type == 'reactions':
                    ranking_matrix = drg.drgep_ranked([case], mec, reduction_type, integrity=integrity)
                elif reduction_type == 'qss':
                    if qss_method == 'AOI':
                        ranking_matrix = qss.qss_ranked_species_AOI([case], mec)
                    elif qss_method == 'LOI':
                        ranking_matrix = qss.qss_ranked_species_LOI([case], mec)
                    else:
                        raise NameError('Wrong value for parameter qss_method !'
                                        'Possible options are : "LOI" or "AOI".')
                elif reduction_type == 'lumping':
                    ranking_matrix = lumping.find_isomers(mec)
                else:
                    raise ValueError('Reduction type not good !')

                # Reordering the dictionary to rank
                ranking_matrix_sorted = sorted(ranking_matrix.items(), key=operator.itemgetter(1), reverse=True)
                ranking_matrix_final = collections.OrderedDict(ranking_matrix_sorted)
                # Construction table
                index = 0
                for id_spec, (spec, coefficient) in enumerate(ranking_matrix_final.items()):
                    index = id_spec - nb_line_max * (id_spec // nb_line_max)
                    if reduction_type == 'lumping':
                        table_data[index].append(str(id_spec + 1) + ' ' + str(coefficient))
                        color_data[index].append(plt.cm.coolwarm(0))
                    else:
                        table_data[index].append(str(id_spec + 1) + ' ' + str(spec) + ' - ' +
                                                 '{:.2e}'.format(coefficient))
                        color_data[index].append(plt.cm.coolwarm(1 / (-np.log(coefficient) + 1)))

                # Reset the maximum number of line for 'lumping' reduction_type
                if reduction_type == 'lumping':
                    if nb_line == 'auto':
                        if nb_lump < len(ranking_matrix_final):
                            nb_lump = len(ranking_matrix_final)

                for id_spec in range(index + 1, nb_line_max, 1):
                    table_data[id_spec].append(' ')
                    color_data[id_spec].append(plt.cm.coolwarm(0.5))

        else:
            # Computing is different if the reduction type is different
            # DRGEP for reaction and species, finding isomers for lumping
            # and LOI for qss.
            if reduction_type == 'species' or reduction_type == 'reactions':
                ranking_matrix = drg.drgep_ranked(case_list, mec, reduction_type, integrity=integrity)
            elif reduction_type == 'qss':
                if qss_method == 'AOI':
                    ranking_matrix = qss.qss_ranked_species_AOI(case_list, mec)
                elif qss_method == 'LOI':
                    ranking_matrix = qss.qss_ranked_species_LOI(case_list, mec)
                else:
                    logger.error('Wrong value for parameter qss_method !')
                    logger.error('Possible options are : "LOI" or "AOI".')
                    exit()
            elif reduction_type == 'lumping':
                ranking_matrix = lumping.find_isomers(mec)
            else:
                ranking_matrix = 0
                logger.error('Reduction type not good !')
            # Reordering the dictionary to rank
            ranking_matrix_sorted = sorted(ranking_matrix.items(), key=operator.itemgetter(1), reverse=True)
            ranking_matrix_final = collections.OrderedDict(ranking_matrix_sorted)
            # Construction table
            index = 0
            for id_spec, (spec, coefficient) in enumerate(ranking_matrix_final.items()):
                index = id_spec - nb_line_max * (id_spec // nb_line_max)
                if reduction_type == 'lumping':
                    table_data[index].append(str(id_spec + 1) + ' ' + str(coefficient))
                    color_data[index].append(plt.cm.coolwarm(0))
                else:
                    table_data[index].append(str(id_spec + 1) + ' ' + str(spec) + ' - ' + '{:.2e}'.format(coefficient))
                    color_data[index].append(plt.cm.coolwarm(1 / (-np.log(coefficient) + 1)))

            # Reset the maximum number of line for 'lumping' reduction_type
            if reduction_type == 'lumping':
                if nb_line == 'auto':
                    if nb_lump < len(ranking_matrix_final):
                        nb_lump = len(ranking_matrix_final)

            for id_spec in range(index + 1, nb_line_max, 1):
                table_data[id_spec].append(' ')
                color_data[id_spec].append(plt.cm.coolwarm(0.5))

    # Return html table for different results
    if reduction_type == 'species':
        show_table(table_data, label_data, directory + '/drgep_species.html')
    elif reduction_type == 'reactions':
        show_table(table_data, label_data, directory + '/drgep_reactions.html')
    elif reduction_type == 'qss':
        show_table(table_data, label_data, directory + '/qss.html')
    elif reduction_type == 'lumping':
        # reshape the size of the table to avoid empty cells (TODO : same with color_data)
        if nb_line == 'auto':
            nb_rm = len(table_data) - nb_lump
            for i in range(nb_rm):
                table_data.remove(table_data[nb_lump])
        show_table(table_data, label_data, directory + '/lumping.html')


def pyrolysis_matrix(case, mechanism):
    """Compute Direct interaction coefficients (DIC) for a given sample

    :param case: Case object
    :param mechanism: Mechanism object

    :return: number_of_species * number_of_species matrix of coefficients for type species
             number_of_species * number_of_reactions matrix of coefficients for type reactions

    Created: 14/11/2017 [PP]
    Last Modified: 03/04/2018 [QC]
    """
    data = case.extract_profile(mechanism)
    names = case.names_dictionary(mechanism)

    sdb = samples_database([case], mechanism)

    maximum_values_dict = {}
    for spec in mechanism.ctmech.species_names:
        value = np.max(data[:, names[spec]])
        if value > 1e-20:
            maximum_values_dict[spec] = np.max(data[:, names[spec]])

    coupled_matrix = np.zeros([mechanism.ns, mechanism.ns])
    for sample in sdb:

        # Parameters
        my_ns = mechanism.ns
        net = mechanism.network
        ct_mechanism = mechanism.ctmech
        nup = net.nup
        nur = net.nur

        # Get reaction rates for each reaction
        ct_mechanism.TPY = float(sample.T), float(sample.P), [float(myY) for myY in sample.Y]

        omega = ct_mechanism.forward_rates_of_progress

        # Get production and consumption rates for species
        pa = ct_mechanism.creation_rates
        ca = ct_mechanism.destruction_rates

        # Workspace
        DIC_spec = np.zeros([my_ns, my_ns], 'd')

        composition = [spec.composition for spec in ct_mechanism.species()]

        # Evaluate DIC(i,j)
        for i in range(my_ns):
            # reactions containing species i
            bool_i = net.indr[net.inds == i]
            ind_j = net.indsj[net.indsi == i]

            for j in ind_j:
                if 'C' in composition[i] and 'C' in composition[j] \
                        and 'H' in composition[i] and 'H' in composition[j]:
                    if composition[j]['C'] <= composition[i]['C'] \
                            and composition[j]['H'] <= composition[i]['H']:
                        if i == j or composition[i] == composition[j]:
                            DIC_spec[i, j] = 0
                        else:
                            bool_j = net.indr[net.inds == j]  # reactions containing species j
                            ind_k = np.intersect1d(bool_i, bool_j)  # reactions containing species i and j

                            # Compute the DIC
                            for k in ind_k:
                                DIC_spec[i, j] += (nup[i, k] - nur[i, k]) * omega[k]

                            # Normalize
                            DIC_spec[i, j] = abs(DIC_spec[i, j]) / max(max(pa[i], ca[i]), 1e-60)

        coupled_matrix = np.fmax(coupled_matrix, DIC_spec)

    return coupled_matrix


def path_matrix(case, mechanism, plot=False):
    """ Compute the path matrix for a given mechanism and case

    :param case: case object
    :param mechanism: mechanism object
    :param plot: boolean to plot or not the path matrix

    :return reduced_matrix, reduced_species_names: double array of coefficient
    and list of strings

    Created: 03/04/2018 [QC]
    """

    data = case.extract_profile(mechanism)
    names = case.names_dictionary(mechanism)

    maximum_values_dict = {}
    for spec in mechanism.ctmech.species_names:
        value = np.max(data[:, names[spec]])
        if value > 1e-20:
            maximum_values_dict[spec] = np.max(data[:, names[spec]])

    species_names = mechanism.ctmech.species_names

    coupled_matrix = pyrolysis_matrix(case, mechanism)

    full_coupled_matrix = coupled_matrix.copy()

    # Discarding too law values
    for i in range(len(species_names)):
        for j in range(len(species_names)):
            if species_names[i] not in maximum_values_dict:
                full_coupled_matrix[i][j] = 0
            if species_names[j] not in maximum_values_dict:
                full_coupled_matrix[j][i] = 0
            # if full_coupled_matrix[i][j] < 1e-10:
            #     full_coupled_matrix[i][j] = np.nan

    import scipy.sparse as sparse
    full_coupled_matrix_sparse = sparse.csr_matrix(full_coupled_matrix)
    array_of_reorder = sparse.csgraph.reverse_cuthill_mckee(full_coupled_matrix_sparse)
    full_coupled_matrix = [[full_coupled_matrix[i][j] for j in array_of_reorder] for i in array_of_reorder]

    full_coupled_matrix = np.array(full_coupled_matrix)

    index_of_zeros = 0
    for i in range(len(species_names)):
        if sum(full_coupled_matrix[i, :]) != 0 or sum(full_coupled_matrix[:, i]) != 0:
            if i > index_of_zeros:
                index_of_zeros = i

    reduced_matrix = full_coupled_matrix[:index_of_zeros, :index_of_zeros]
    reduced_species_names = [species_names[i] for i in array_of_reorder[:index_of_zeros]]

    plt.figure(num=None, figsize=(10, 7.5), dpi=100, facecolor='w', edgecolor='k')
    plt.pcolor(reduced_matrix, cmap='RdBu')
    x = list(range(len(reduced_species_names)))
    labels = reduced_species_names
    plt.xticks(x, labels, rotation='vertical')
    plt.yticks(x, labels)
    plt.tick_params(labelsize=10)
    plt.colorbar()
    if plot:
        plt.show()
    else:
        plt.savefig('path_matrix', dpi=100)

    return reduced_matrix, reduced_species_names


def coupling_matrix_graph(mechanism):
    """ Return the coupling matrix for plotting graph

    :param mechanism: mechanism object

    :return:

    Created: 03/04/2018 [QC]
    """

    species_names = mechanism.ctmech.species_names

    coupled_matrix = nk.Network(mechanism).deltaSS

    full_coupled_matrix = coupled_matrix.copy()
    import scipy.sparse as sparse
    full_coupled_matrix_sparse = sparse.csr_matrix(full_coupled_matrix)
    array_of_reorder = sparse.csgraph.reverse_cuthill_mckee(full_coupled_matrix_sparse)
    full_coupled_matrix = [[full_coupled_matrix[i][j] for j in array_of_reorder] for i in array_of_reorder]
    species_names = [species_names[i] for i in array_of_reorder]

    for i in range(len(species_names)):
        for j in range(len(species_names)):
            if full_coupled_matrix[i][j] == 1 and full_coupled_matrix[j][i] == 1:
                full_coupled_matrix[i][j] = 2
                full_coupled_matrix[j][i] = 2

    plt.figure(num=None, figsize=(10, 7.5), dpi=100, facecolor='w', edgecolor='k')
    plt.pcolor(full_coupled_matrix, cmap='RdBu')
    x = list(range(len(species_names)))
    labels = species_names
    plt.xticks(x, labels, rotation='vertical')
    plt.yticks(x, labels)
    plt.tick_params(labelsize=10)
    plt.colorbar()
    plt.savefig('coupled_matrix', dpi=100)


def pyrolysis_products(case, mechanism, fuel_species, plot=False):
    """ Return the pyrolysis products of a given fuel species

    :param case: case object
    :param mechanism: mechanism object
    :param fuel_species: string representing the fuel species
    :param plot: boolean to plot the path matrix or not

    :return:

    Created: 03/04/2018 [QC]
    """

    matrix, names = path_matrix(case, mechanism, plot=plot)

    if type(fuel_species) != list:
        fuel_species = [fuel_species]

    products_by_level = {0: fuel_species}
    new_products_of_level = {0: fuel_species}
    all_species = []

    level = 1
    while new_products_of_level[level - 1]:
        products_by_level[level] = []
        new_products_of_level[level] = []

        all_species += products_by_level[level - 1]
        all_species = list(set(all_species))

        for spec in products_by_level[level - 1]:
            index_spec = names.index(spec)
            products_of_species = [spec for index, spec in enumerate(names)
                                   if matrix[index_spec, index] > 0]

            products_by_level[level] += products_of_species

        products_by_level[level] = list(set(products_by_level[level]))

        for spec in products_by_level[level]:
            if spec not in all_species:
                new_products_of_level[level].append(spec)

        level += 1

    logger.info('\n### Species produced at each level from fuel species ###')
    for level in products_by_level:
        logger.info('Level ' + str(level) + ' products ---> ' + str(products_by_level[level]) + '='
                    + str(len(products_by_level[level])))

    logger.info('\n### New species appearing at each level from fuel species ###')
    for level in products_by_level:
        logger.info('Level ' + str(level) + ' new products ---> ' + str(new_products_of_level[level])
                    + '=' + str(len(new_products_of_level[level])))
        if len(new_products_of_level[level]) == 0:
            break

    logger.info('\nTotal of species involved = ' + str(len(all_species)))


def find_decomposition(mechanism, species_of_interest):
    """ Finds the reactions involved in the decomposition of a species

    :param mechanism: mechanism object
    :param species_of_interest: string representing the species of interest

    :return: decomposition_reactions, decomposition_species

    Created: 30/08/2019 [QC]
    """

    ct_mechanism = mechanism.ctmech

    species = ct_mechanism.species()
    species_names = ct_mechanism.species_names
    reactions = ct_mechanism.reactions()

    decomposition_reactions = []

    decomposition_species = []

    added_reaction = True
    added_species = True

    while added_reaction and added_species:

        for reaction in reactions:
            reactants = reaction.reactants
            products = reaction.products

            for reactant in reactants:
                if 'C' in species[species_names.index(reactant)].composition:
                    if reactant in decomposition_species or reactant in species_of_interest:
                        decomposition_reactions.append(reaction)
                        added_reaction = True
                    else:
                        added_reaction = False

                    for product in products:
                        if 'C' in species[species_names.index(product)].composition:
                            if species[species_names.index(reactant)].composition['C'] >= \
                                    species[species_names.index(product)].composition['C'] > 4 \
                                    and product not in decomposition_species:
                                decomposition_species.append(product)
                                added_species = True
                            else:
                                added_species = False

    return decomposition_reactions, decomposition_species


def hr_analysis(case, method='species', quantity='HR', integrated=False, movie=True,
                directory=None, size='FullHD', threshold=20, fmt='jpg', plot=False,
                tstart=None, tend=None):
    """ Function to draw a bar plot on which species or reaction contributes the
    most to the creation of heat release for a given mechanism and a given
    case.

    :param case: case object
    :param method: method used, either species or reactions
    :param quantity: target quantity for the analysis ("HR" or "Qj")
    :param integrated: keyword to consider the classical sample database or the
    integrated one
    :param movie: create a movie or not (boolean)
    :param directory: choose the directory where the images and movies are
    created (string)
    :param size: control the size of the movie (string)
    :param threshold: control the number of variable that are going to be
    displayed on the bar plot
    :param fmt: output format
    :param plot: boolean to enable plotting at the end of the function or not
    :param tstart: first position/time to start the analysis(float) 
    :param tend: last position/time to end the analysis(float)

    :return hr_element: panda datatable of hr
    :return hr_element_int: panda datatable of hr integrated in the interval [tsart ; tend]

    Created: 29/08/2019 [JW]
    Last modification: 25/06/2020 [AP]
        - It is possible to integrate the analysis over multiple samples
        - Only the samples in the interval define by [tsart ; tend] will be analysed
        - Add Qj information to the HR analysis to know the reaction way (backward or forward) 
        - Add digits to the filename

    """

    logger.info("Performing analysis on HR")
    logger.info("--------------------------")

    # Create the folder where the image + video are gonna be put in
    if not directory:
        directory = analysis_dir + '/analysis_' + quantity + '_folder/'
    database.create_dir(directory, overwrite=False)

    images = []
    importance = 0.4
    # Check size of the plot
    if size == 'QHD':
        size_vector = [2560, 1440]
    elif size == 'FullHD':
        size_vector = [1920, 1080]
    elif size == 'HD':
        size_vector = [1280, 720]
    elif size == 'qHD':
        size_vector = [960, 540]
    elif size == 'YouTube144p':
        size_vector = [256, 144]
    elif isinstance(size, list):
        size_vector = [size[0], size[1]]
    else:
        logger.error('ERROR ! You should specify either the size of the vector or the keywords '
                     'QHD, FullHD, HD, qHD, Youtube144p.')
        sys.exit()

    # Check quantity
    if method == "species" and quantity == "Qj":
        logger.error('\nERROR ! The Qj analysis is not available for species.')
        logger.error('==> The quantity parameter has been set to HR.\n')

    # Check integrated intervals
    if not tstart:
        tstart = 0.0
    if not tend:
        tend = 10000.0

    # Create sample database on cases
    logger.info("=> Creating sampling database")
    logger.info("#############################")
    sdb = samples_database([case], case.mechanism, integrated)

    # Get solution file, reaction names and species names
    ct_mechanism = case.mechanism.ctmech
    if case.mechanism.f90:
        species = ct_mechanism.skeletal.species_names
        reactions = ct_mechanism.skeletal.reaction_equations()
    else:
        species = ct_mechanism.species_names
        reactions = ct_mechanism.reaction_equations()

    if method == 'species':
        element = species
    elif method == 'reactions':
        element = reactions
    else:
        raise NameError('The method is not right, it should be either species or reactions.')

    pd.options.display.max_rows = len(element)
    pd.options.display.max_columns = len(sdb)

    hr_element = pd.DataFrame(data=[], index=element)
    hr_element_int = pd.DataFrame(data=[], index=element)
    first_int = True

    for id_sample, my_sample in enumerate(sdb):
        if tstart <= my_sample.grid <= tend:
            if method == 'reactions' and quantity == "HR":
                # list with name of reactions and reaction rate of progress (reset for each sample)
                reactions_QJ = []

            logger.info('Computing sample ' + str(id_sample) + ' for grid point ' + str(my_sample.grid))
            hr_element[str(id_sample)] = ""

            # Get HR or Qj
            if quantity == "HR":
                hr_element_sample = globals()[method + "_heat_release"](case.mechanism, my_sample)
                ct_mechanism.TPY = float(my_sample.T), float(my_sample.P), [float(myY) for myY in my_sample.Y]
                QJ = ct_mechanism.net_rates_of_progress

            elif quantity == "Qj":
                if case.mechanism.f90:
                    hr_element_sample = tools.compute_rates(case.mechanism, my_sample.T, my_sample.P,
                                                            [float(myY) for myY in my_sample.Y])
                else:
                    ct_mechanism.TPY = float(my_sample.T), float(my_sample.P), [float(myY) for myY in my_sample.Y]
                    hr_element_sample = ct_mechanism.net_rates_of_progress
            else:
                raise NameError('Quantity does not exist. It should be either HR or Qj.')

            # Compute integrated variables
            if first_int:
                hr_element_sample_int = np.zeros(len(hr_element_sample))
                if method == 'reactions' and quantity == "HR":
                    QJ_int = np.zeros(len(QJ))
                old_grid = my_sample.grid
                first_int = False
            else:
                for ind_hr, hr in enumerate(hr_element_sample):
                    hr_element_sample_int[ind_hr] += hr * (my_sample.grid - old_grid)
                    if method == 'reactions' and quantity == "HR":
                        QJ_int[ind_hr] += QJ[ind_hr] * (my_sample.grid - old_grid)
                old_grid = my_sample.grid

            # Fill Panda database for hr element 
            for ind_hr, hr in enumerate(hr_element_sample):
                hr_element[str(id_sample)][ind_hr] = hr
                if method == 'reactions' and quantity == "HR":
                    # Add the reaction rate of progress
                    reactions_QJ.append(reactions[ind_hr] + '\n Qj = ' + str(QJ[ind_hr]))

            # Change the index name in the panda dataframe to add reaction rate of progress
            if method == 'reactions' and quantity == "HR":
                hr_element.index = reactions_QJ

            # Construction of the graph
            if threshold > len(element):
                threshold = len(element)

            # Sort values by maximum
            indices = hr_element[str(id_sample)].abs().sort_values(ascending=False).index

            if plot:
                if quantity == "HR":
                    hr_element[str(id_sample)].loc[indices].iloc[0:threshold].plot.barh(
                        title="Contribution to the heat release by the " + method +
                              ": \n hr=" + '{0:.2e}'.format(sum(hr_element_sample)), figsize=(10, 10))
                    plt.xlabel(r'Heat released by ' + method + ' : $w_{kj}$')

                elif quantity == "Qj":
                    hr_element[str(id_sample)].loc[indices].iloc[0:threshold].plot.barh(
                        title="Maximum reaction rates of progress \n T = " + '{:.1f}'.format(my_sample.T) + "K",
                        figsize=(10, 10))
                    plt.xlabel(r'Reaction rates of progress: $Q_{j}$')
                plt.gca().invert_yaxis()
                plt.tight_layout()
                plt.savefig(directory + '/graph.' + fmt)
                plt.close()

                # Create aside graph showing the physics
                logger.info("=> Drawing graph for case")
                ready_made_plot(directory, 'plot1.' + fmt, integrated, ct_mechanism, sdb, my_sample, importance,
                                size_vector, zoom=False)

                # Create aside zoom graph showing the physics
                logger.info("=> Drawing graph for case (near view)")
                ready_made_plot(directory, 'plot2.' + fmt, integrated, ct_mechanism, sdb, my_sample, importance,
                                size_vector, zoom=True)

                # Merge the two plot graphs (physic of the flame + zoom)
                logger.info("=> Merging graphs")
                x1, x2 = int(size_vector[0] * importance), int(size_vector[0] * importance)
                y1, y2 = int(size_vector[1] / 2), int(size_vector[1] / 2)
                img_merge = merge_vertical('plot1.' + fmt, x1, y1, 'plot2.' + fmt, x2, y2, directory)
                img_merge.save(directory + '/plot.' + fmt)

                # Create merge of the two graphs (previous merge + digraph)
                x1, x2 = int(size_vector[0] * (1 - importance)), int(size_vector[0] * importance)
                y1, y2 = int(size_vector[1]), int(size_vector[1])
                img_merge = merge_horizontal('graph.' + fmt, x1, y1, 'plot.' + fmt, x2, y2, directory)
                img_merge.save(directory + '/Prod_' + method + '_path_' + '{:.10f}'.format(my_sample.grid) + '.' + fmt)
                images.append('Prod_' + method + '_path_' + '{:.10f}'.format(my_sample.grid) + '.' + fmt)

    # Construction of the database for integrated values
    hr_element_int["0"] = ""
    reactions_QJ = []
    for ind_hr, hr in enumerate(hr_element_sample_int):
        hr_element_int["0"][ind_hr] = hr

        if method == 'reactions' and quantity == "HR":
            # Add the reaction rate of progress
            reactions_QJ.append(reactions[ind_hr] + '\n' + r" $\overline{Q_j}$ = "
                                + str(QJ_int[ind_hr] / (old_grid - tstart)))

    # Change the index name in the panda dataframe to add reaction rate of progress
    if method == 'reactions' and quantity == "HR":
        hr_element_int.index = reactions_QJ

    # Construction of the integrated histogram
    # Sort values by maximum
    indices = hr_element_int["0"].abs().sort_values(ascending=False).index
    if plot:
        if quantity == "HR":
            hr_element_int["0"].loc[indices].iloc[0:threshold].plot.barh(
                title="Integrated contribution to the heat release " +
                      "\n t/x = [{tstart:.2e} ; {tend:.2e}]".format(tstart=tstart, tend=tend), figsize=(10, 10))
            plt.xlabel(r'Integrated heat released')

        elif quantity == "Qj":
            hr_element_int["0"].loc[indices].iloc[0:threshold].plot.barh(
                title="Maximum integrated reaction rates of progress \n t/x = [{tstart:.2e} ; {tend:.2e}]".format(
                    tstart=tstart, tend=tend),
                figsize=(10, 10))
            plt.xlabel(r'Integrated reaction rates of progress')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(directory + '/graph_integrated.' + fmt)
        plt.close()

    if plot:
        os.remove(directory + '/graph.' + fmt)
        os.remove(directory + '/plot.' + fmt)

    # Send message to console
    indices = hr_element.abs().max(axis=1).sort_values(ascending=False).index
    logger.info(str(hr_element.loc[indices].head(len(element))))
    logger.info("=> Integrated results")
    indices = hr_element_int["0"].abs().sort_values(ascending=False).index
    logger.info(str(hr_element_int.loc[indices].head(len(element))))
    # Create the movie if asked
    if plot and movie:
        logger.info("=> Creating movie")
        from cv2 import VideoWriter, VideoWriter_fourcc, imread, resize

        # Chose the folder to save the movie
        image_folder = directory
        video_name = directory + 'Reaction_prod.avi'

        # Write the movie
        fourcc = VideoWriter_fourcc(*"XVID")
        video = VideoWriter(video_name, fourcc, float(1), (size_vector[0], size_vector[1]))
        for image in images:
            video.write(imread(os.path.join(image_folder, image)))

        # Save the movie
        video.release()

        logger.info("=> Movie created")

    logger.info("=> Done")

    return hr_element, hr_element_int


def sensitivity_analysis(case, reduction_type, method='perturb', perturbation=5e-4, threshold=0, drgep=False,
                         omega=False, overwrite=False, restore='no', restore_data=True, panda=False, plot=False,
                         directory=None):
    """ Function to analyse the sensitivity analysis of a mechanism.
    Adapted from Cantera open source website.
    Advised value for perturbation is 5e-4 but you may increase it (if the
    simulation crashes) or decrease it (if there are no changes).

    :param case: list of case object
    :param reduction_type: string of reduction type
    :param method: string of method type
    :param perturbation: perturbs the reaction rate coefficients (double)
    :param threshold: value of the sensor to consider changing value of DRG or not
    :param drgep: enable the ranking according to the DRGEP
    :param omega: parameter to distinguish positive and negative effects
    :param overwrite: boolean to decide whether cases should be overwritten or not
    :param restore_data: boolean to enable the code to restore data from previous sensibility analysis
    :param panda: boolean to decide whether to sort it as panda database or dictionary
    :param plot: plot the sensitivities if wanted
    :param directory: string to specify a user-defined directory

    :return sa : import sensitive species or reactions

    Created : 30/08/2019 [JW]
    """

    # Create the folder where the sensibility restore file is gonna be put in
    if not directory:
        directory = analysis_dir + '/analysis_sensitivity/'
    database.create_dir(directory, overwrite=False)

    # Case mechanism
    # For ARC mechanism, be sure that the mechanism is well recompiled to avoid numerical errors.
    mechanism = case.mechanism
    if case.mechanism.f90:
        case.mechanism.reset()
        f90_filename = mechanism.mechdir + "/" + str(mechanism.name) + '.f90'
        custom_kinetics.print_fortran(case.mechanism,
                                      f90_filename=f90_filename,
                                      species_qss_names=[], use='Cantera', routine_name='customkinetics',
                                      semi_implicit=False, rrate=True, force=True, exponential=False,
                                      maintain_given_order=tools.convert_to_valid_fortran_var
                                      (case.mechanism.species_qss_names[::-1]))
        custom_kinetics.compile_fortran(f90_filename, force=False, fortran_format='f90', mechanism_lib=None,
                                        output=False, remote_install=True)
        case.mechanism.f90 = f90_filename
        case.mechanism.ctmech = ct.Solution(case.mechanism.path)
    gas = mechanism.ctmech

    # Set species and reactions lists
    if case.mechanism.f90:
        species = mechanism.skeletal.ctmech.species_names
        equations = mechanism.skeletal.ctmech.reactions()
    else:
        species = gas.species_names
        equations = gas.reactions()

    logger.info("Mechanism is " + str(mechanism.name))
    logger.info("---------------------")

    # Gather length of mechanism list, case list and error list to generate panda dataframe
    case_err_list = []
    case_id_list = []
    case_mechanism_list = []
    for err in list(case.error_dict.keys()):
        case_err_list.append(err)
        case_id_list.append(case.myid)
        case_mechanism_list.append(case.mechanism.name)
    columns = pd.MultiIndex.from_arrays([case_mechanism_list, case_id_list, case_err_list])

    # Separate cases where reduction_type is species and reactions
    if reduction_type in ['species', 'S']:
        n_loop = len(species)
        sensitivities = pd.DataFrame(data=[], index=range(len(species)), columns=columns)
    elif reduction_type in ['reactions', 'R']:
        n_loop = len(equations)
        sensitivities = pd.DataFrame(data=[], index=range(len(equations)), columns=columns)
    else:
        raise NameError("Problem with the reduction type !")

    # Set showing parameters (otherwise matrix does not show completely to the screen)
    pd.options.display.max_rows = max(len(species), len(equations))
    pd.options.display.max_columns = len(case_err_list) + 1

    # Threshold above 1 means a certain number of entities displayed and below 1 a threshold on the sensitivity
    # coefficient. Threshold equal 1 means nothing is done.
    if threshold > 1:
        threshold_number = int(threshold)
        threshold_value = 0.0
        logger.info("Performing analysis on the sensitivity with a perturbation of " + str(perturbation)
                    + " and for the biggest " + str(threshold_number) + ".")
    elif threshold < 1:
        threshold_value = threshold
        threshold_number = n_loop
        logger.info("Performing analysis on the sensitivity with a perturbation of " + str(perturbation)
                    + " and a threshold of " + str(threshold_value))
        logger.info("--------------------------------------")
    else:
        threshold_value = 0.0
        threshold_number = n_loop
        logger.info("Performing analysis on the sensitivity with a perturbation of " + str(perturbation)
                    + " and no threshold.")
        logger.info("--------------------------------------")

    # Run case of reference
    logger.info("\n ### Run reference cases")
    case.run(mechanism, loglevel=1, overwrite=overwrite, restore=restore)
    if not case.success:
        logger.error("\n/!\\ ERROR /!\\ The reference case for sensibilities failed,"
                     "then sensibilities cannot be computed")
        logger.error("==> EXIT")
    if case.mechanism.f90:
        case.mechanism.reset()

    # Compute DRGEP and rank reactions
    if drgep:
        if reduction_type in ['species', 'S']:
            logger.info("\n ### Compute DRGEP for species")
            drgep_dict = drg.drgep_ranked([case], case.mechanism, reduction_type='S')
        elif reduction_type in ['reactions', 'R']:
            logger.info("\n ### Compute DRGEP for reactions")
            drgep_dict = drg.drgep_ranked([case], case.mechanism, reduction_type='R')
        else:
            raise NameError('Reduction type is not right !')

        ranked = sorted(zip(drgep_dict.keys(), drgep_dict.values()), reverse=True, key=lambda x: x[1])
        drgep_elements = ranked[0:threshold_number]
        for element in ranked:
            if element[1] < threshold_value:
                drgep_elements.remove(element)

    # Computation of the rate (net_production_rate if species and progress_rate if reactions) at maximum HR.
    # Can be useful if the user wants to correlate sign of sensitivity and sign of rate.
    if omega:
        data_dict = case.data_dict(mechanism)
        ind_max_hr = np.argmax(data_dict['HR'][:])
        gas_max_hr = mechanism.ctmech
        gas_max_hr.TPY = data_dict['T'][ind_max_hr], data_dict['P'][ind_max_hr], data_dict['Y'][ind_max_hr]
        if reduction_type in ['species', 'S']:
            omega_matrix = gas_max_hr.net_production_rates
        elif reduction_type in ['reactions', 'R']:
            omega_matrix = gas_max_hr.net_rates_of_progress

    # Case filename and restore data if allowed and if any
    case_filename = case.myid + 'method' + str(method) + 'perturbation' + str(perturbation)
    restore_filename = directory + case.mechanism.name + case.myid + '_type_' + reduction_type \
                       + '_method_' + str(method) + '_perturbation_' + str(perturbation)
    file_restore = False
    if os.path.isfile(restore_filename) and not overwrite and restore_data:
        sensitivities = pd.read_pickle(restore_filename)
        file_restore = True

    # Loop over the interesting quantity
    for m in range(n_loop):
        if file_restore:
            break
        # Create another mechanism to change reaction coefficients
        # Method 1 is the perturbation of the pre-Arrhenius constant.
        if method == 'perturb':
            if reduction_type in ['species', 'S']:
                # if drgep, treats only the cases that have been selected previously
                if drgep and species[m] not in [elt[0] for elt in drgep_elements]:
                    continue

                # copy the case and reinitialize the data stored.
                case_change = case.copy(case_filename + 'species' + str(species[m]) + '_' + str(m))
                case_change.stored_data = {}

                # If f90, print the new f90 and compile.
                if case.mechanism.f90:
                    case.mechanism.reset()
                    case_change.mechanism.reset()
                    for id_eq, eq in enumerate(equations):
                        if species[m] in eq.reactants.keys() or species[m] in eq.products.keys():
                            case_change.mechanism.skeletal.ctmech.set_multiplier(1 + perturbation, id_eq)
                    case_change.mechanism.parent = case_change.mechanism.skeletal
                    f90_filename = mechanism.mechdir + "/" + str(mechanism.name) + '_' + species[m] + '_' + str(
                        m) + '.f90'
                    custom_kinetics.print_fortran(case_change.mechanism, f90_filename=f90_filename,
                                                  species_qss_names=[], use='Cantera', routine_name='customkinetics',
                                                  semi_implicit=False, rrate=True, force=True, exponential=False,
                                                  maintain_given_order=tools.convert_to_valid_fortran_var
                                                  (case_change.mechanism.species_qss_names[::-1]))
                    custom_kinetics.compile_fortran(f90_filename, force=False, fortran_format='f90', mechanism_lib=None,
                                                    output=False, remote_install=True)
                    case_change.mechanism.f90 = f90_filename
                    case_change.mechanism.ctmech = ct.Solution(case.mechanism.path)
                else:
                    for id_eq, eq in enumerate(equations):
                        if species[m] in eq.reactants.keys() or species[m] in eq.products.keys():
                            case_change.mechanism.ctmech.set_multiplier(1 + perturbation, id_eq)

                logger.info("------------------")
                logger.info("Perturbing species " + str(species[m]))
            elif reduction_type in ['reactions', 'R']:
                # if drgep, treats only the cases that have been selected previously
                if drgep and equations[m].equation not in [elt[0].equation for elt in drgep_elements]:
                    continue
                # copy the case and reinitialize the data stored.
                case_change = case.copy(case_filename + 'reactions' + str(equations[m]) + '_' + str(m))
                case_change.stored_data = {}

                # If f90, print the new f90 and compile.
                if case.mechanism.f90:
                    case.mechanism.reset()
                    case_change.mechanism.reset()
                    case_change.mechanism.skeletal.ctmech.set_multiplier(1 + perturbation, m)
                    case_change.mechanism.parent = case_change.mechanism.skeletal
                    f90_filename = mechanism.mechdir + "/" + str(mechanism.name) + '_' + str(equations[m]) + '_' + str(
                        m) \
                                   + '.f90'
                    custom_kinetics.print_fortran(case_change.mechanism, f90_filename=f90_filename,
                                                  species_qss_names=[], use='Cantera', routine_name='customkinetics',
                                                  semi_implicit=False, rrate=True, force=True, exponential=False,
                                                  maintain_given_order=tools.convert_to_valid_fortran_var
                                                  (case_change.mechanism.species_qss_names[::-1]))
                    custom_kinetics.compile_fortran(f90_filename, force=False, fortran_format='f90', mechanism_lib=None,
                                                    output=False, remote_install=True)
                    case_change.mechanism.f90 = f90_filename
                    case_change.mechanism.ctmech = ct.Solution(case.mechanism.path)
                else:
                    case_change.mechanism.ctmech.set_multiplier(1 + perturbation, m)  # perturb reaction m

                logger.info("--------------------")
                logger.info("Perturbing reaction " + str(equations[m]))
            else:
                raise NameError('Reduction type is not right !')

        # Method 2 is the creation of another mechanism with deleted species or reactions
        elif method == 'delete':
            if reduction_type in ['species', 'S']:
                # copy the case and reinitialize the data stored.
                case_change = case.copy(case_filename + 'species' + str(species[m]))
                case_change.stored_data = {}

                # change the mechanism with the deleted species or reaction
                ctmech_change = custom_kinetics.cantera_reduced_solution(mechanism.ctmech,
                                                                         species_to_discard=species[m])
                reduced_how = 'rmS'
                case_change.mechanism = mechanisms.Mechanism(parent=mechanism, cti=ctmech_change, how=reduced_how)

            elif reduction_type in ['reactions', 'R']:
                # copy the case and reinitialize the data stored.
                case_change = case.copy(case_filename + 'reactions' + str(equations[m]))
                case_change.stored_data = {}

                # change the mechanism with the deleted species or reaction
                ctmech_change = custom_kinetics.cantera_reduced_solution(mechanism.ctmech,
                                                                         reactions_to_discard=[equations[m]])
                reduced_how = 'rmR'
                case_change.mechanism = mechanisms.Mechanism(parent=mechanism, cti=ctmech_change, how=reduced_how)
            else:
                raise NameError('Reduction type is not right !')

        else:
            raise NameError('The method can be either perturb or delete.')

        # Compute the new case
        try:
            case_change.run(case_change.mechanism, loglevel=0, overwrite=overwrite, restore=restore)
        except RuntimeError:
            logger.warning('The case has not converged.')
            pass

        if case_change.mechanism.f90:
            case.mechanism.reset()
            case_change.mechanism.reset()

        # Compute sensitivity coefficient (X_change - X)/(X*dk) with X the target variable of the case
        relative_error = error.case_error(case, mechanism, case=case_change, mechanism=case_change.mechanism,
                                          error_factor=1)
        for id_err, err in enumerate(relative_error):
            if reduction_type in ['species', 'S']:
                if omega and omega_matrix[m] < 0:
                    sensitivities.loc[m, mechanism.name][case.myid][list(case.error_dict.keys())[id_err]] \
                        = - err / perturbation
                elif omega and omega_matrix[m] < 0:
                    sensitivities.loc[m, mechanism.name][case.myid][list(case.error_dict.keys())[id_err]] \
                        = err / perturbation
                else:
                    sensitivities.loc[m, mechanism.name][case.myid][list(case.error_dict.keys())[id_err]] \
                        = abs(err / perturbation)

            elif reduction_type in ['reactions', 'R']:
                if omega and omega_matrix[m] < 0:
                    sensitivities.loc[m, mechanism.name][case.myid][list(case.error_dict.keys())[id_err]] \
                        = - err / perturbation
                elif omega and omega_matrix[m] < 0:
                    sensitivities.loc[m, mechanism.name][case.myid][list(case.error_dict.keys())[id_err]] \
                        = err / perturbation
                else:
                    sensitivities.loc[m, mechanism.name][case.myid][list(case.error_dict.keys())[id_err]] \
                        = abs(err / perturbation)

        if omega and omega_matrix[m] < 0:
            logger.info("!!!WARNING!!!")
            logger.info("This reaction is reverse, the sign of the sensitivity has been changed ")

    # Rank coefficients
    indices = sensitivities[case.mechanism.name, case.myid].abs().max(axis=1).sort_values(ascending=False).index

    # Display the sensitivity coefficient table
    logger.info("Resume of the table")
    logger.info("-------------------")
    if reduction_type in ['species', 'S']:
        logger.info(str(sensitivities.loc[indices].head(len(species))))
    elif reduction_type in ['reactions', 'R']:
        logger.info(str(sensitivities.loc[indices].head(len(species))))

    # Plot if wanted
    if plot:
        sensitivities[case.mechanism.name, case.myid].loc[indices[0:threshold_number]].dropna() \
            .plot.bar(title="Sensitivities for \n Mechanism : " + str(case.mechanism.name)
                            + " \n Case : " + str(case.myid))

        plt.gca().invert_yaxis()
        plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$')
        plt.tight_layout()
        plt.show()

    # Output on panda format or dictionary format
    sensitivities.to_pickle(restore_filename)
    if panda:
        return sensitivities
    else:
        final_sensitivities = {mechanism.name: {case.myid: {}}}
        for err in list(case.error_dict.keys()):
            final_sensitivities[mechanism.name][case.myid][err] = {}
        for m in range(n_loop):
            for id_err, err in enumerate(list(case.error_dict.keys())):
                if reduction_type in ['species', 'S']:
                    final_sensitivities[case.mechanism.name][case.myid][err][m] = \
                        sensitivities[case.mechanism.name, case.myid, err][m]
                elif reduction_type in ['reactions', 'R']:
                    final_sensitivities[case.mechanism.name][case.myid][err][m] = \
                        sensitivities[case.mechanism.name, case.myid, err][m]
        return final_sensitivities


def uncertainty_analysis(case, mechanism, bounds=None, n_sample=1000, threshold=20, save=False, directory=None):
    """ Function to analyse the uncertainty of a mechanism.

    :param case: case object
    :param mechanism: mechanism object
    :param bounds: set the boundary of the reaction uncertainties (size is
    2*n_reactions)
    :param n_sample: number of sample of the uncertainty analysis
    :param threshold: threshold on the number of elements plotted
    :param save: save the figure or not (boolean)
    :param directory: string to specify the directory to store file

    :return uncertainties: uncertainty matrices withe coefficients

    Created: 29/08/2019 [JW]

    """

    # Create the folder where the sensibility restore file is gonna be put in
    if not directory:
        directory = analysis_dir + '/analysis_uncertainty/'
    database.create_dir(directory, overwrite=False)

    # Useful import SALib
    from SALib.sample import saltelli
    from SALib.analyze import sobol

    logger.info("Performing analysis on the uncertainty")
    logger.info("--------------------------------------")

    # Useful initialisation
    equation_buffer = []
    equation_case = []

    gas = mechanism.ctmech

    pd.options.display.max_rows = gas.n_reactions
    pd.options.display.max_columns = 2

    # Equations are renamed to avoid duplicate errors
    for eq in gas.reaction_equations(range(gas.n_reactions)):
        if eq not in equation_buffer:
            equation_buffer.append(eq)
        else:
            equation_buffer.append(eq + '-duplicate')
    equation_case.append(equation_buffer)

    # Set up boundaries if not set up by user
    if not bounds:
        for n in range(gas.n_reactions):
            bounds.append([0.8, 1.2])

    # Set up problem for sampling and analysis
    problem = {'num_vars': gas.n_reactions, 'names': equation_case, 'bounds': bounds}

    # Sample, evaluate function and analyse variances with Sobol
    param_values = saltelli.sample(problem, n_sample)
    mass_fraction = evaluate(param_values, case, mechanism)
    si_tig = sobol.analyze(problem, mass_fraction[:, 1])
    si_temperature = sobol.analyze(problem, mass_fraction[:, 0])

    # Create uncertainty panda database
    uncertainties = pd.DataFrame(data=[], index=equation_case)
    uncertainties['Tig'] = si_tig['S1'] * 100
    uncertainties['$T_{max}$'] = si_temperature['S1'] * 100

    # For plotting, collect only those steps that are above the threshold
    # Otherwise, the y-axis gets crowded and illegible
    indices = uncertainties.abs().max(axis=1).sort_values(ascending=False).index
    logger.info(str(uncertainties.loc[indices].head(gas.n_reactions)))
    if threshold > gas.n_reactions:
        threshold = gas.n_reactions

    # Create heat map graph
    plt.pcolor(uncertainties.loc[indices].iloc[0:threshold], cmap='RdBu')
    plt.yticks(np.arange(0.5, len(uncertainties.index), 1), uncertainties.index)
    plt.xticks(np.arange(0.5, len(uncertainties.columns), 1), uncertainties.columns)
    plt.tick_params(labelsize=10)
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()

    # Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
    if save is True:
        plt.savefig(directory + 'uncertaintyPlot', dpi=300)

    return uncertainties


def evaluate(values, case, mechanism):
    """ Function to evaluate a parameter with given values.

    :param values: values given by the sampling
    :param case: case object
    :param mechanism: mechanism object

    :return solution: solution vector containing special variable and maximum of
    temperature

    Created: 29/08/2019 [JW]

    """

    # Create solution vector and set up gas object
    solution = np.zeros((values.shape[0], 2))
    gas = mechanism.ctmech

    # Iterate through sampling
    for i, X in enumerate(values):
        # Set the multiplier to the sampling coefficient and run case
        for n in range(gas.n_reactions):
            gas.set_multiplier(n, X[n])
            case.run(mechanism, loglevel=0, overwrite=True, restore='auto')
        # Calculate interesting variable : Tig for 0D, Sl for 1D premixed and
        # integral if HR for diffusion flames and maximum of temperature for
        # all the cases
        if case.reactor_type == 'C0DV' or case.reactor_type == 'C0DP':
            time = case.extract_quantity('time', mechanism=mechanism)
            hr = case.extract_quantity('HR', mechanism=mechanism)
            temperature = case.extract_quantity('T', mechanism=mechanism)
            solution[i, 0] = tools.get_ignition_delay_time(case, mechanism)
            solution[i, 1] = np.max(temperature)
        elif case.reactor_type == 'C1DP':
            solution[i, 0] = case.extract_quantity('Sl', mechanism=mechanism)
            solution[i, 1] = case.extract_quantity('T max', mechanism=mechanism)
        elif case.reactor_type == 'C1CDF' or case.reactor_type == 'C1DF':
            hr = case.extract_quantity('HR', mechanism=mechanism)
            x = case.extract_quantity('Grid', mechanism=mechanism)
            solution[i, 0] = np.trapz(hr, x)
            solution[i, 1] = case.extract_quantity('T max', mechanism=mechanism)
        # reset multipliers
        for n in range(gas.n_reactions):
            gas.set_multiplier(n, 1)

    return solution


def path_analysis(case_list, mechanism, spec, analys='StoS', method='ATOM', integrated=False, movie=True,
                  directory=None, size='FullHD', parse=1, threshold=10 ** -1, atom='C', t_init=None,
                  t_end=None, fmt='jpg', spec2plot=None):
    """ Function to draw fluxes of species either with a Path Flux Analysis
    method or with an atom flux analysis method. It is possible to create
    a species to species graph, a reaction to species graph and a species to
    reaction to species graph.

    :param case_list: list of case object
    :param mechanism: mechanism object
    :param spec: species to consider. Can be a list of species for StoRtoS (string or list of strings)
    :param analys: type of analysis, either StoS (one species related to other species),
    RtoS (one species related to other reactions) or StoRtoS (several species related to several reactions) (string)
    :param method: method used, either PFA (path flux analysis) or ATOM (atom flux analysis) (string)
    :param integrated: keyword to consider the integrated sample database or every point (boolean)
    :param movie: create a movie or not (boolean)
    :param directory: choose the directory where the images and movies are created (string)
    :param size: control the size of the movie (string)
    :param parse: choose sample every parse point
    :param threshold: control the minimum size of the fluxes for which the path is
    showed on the graph
    :param atom: Atom used for the atom flux method (string)
    :param t_init: Initial times or positions for integrated graph (list of floats)
    :param t_end: Final times or positions for integrated graph (list of floats)
    :param fmt: Type of output format (jpg)
    :param spec2plot: Name of the species to display in the graph il "all" is specified
    all the species will be displayed (list of strings)

    :return None

    Created: 29/08/2019 [JW]
    Last modification 29/04/2020 [AP]
      - Add the option to chose the species to display
      - Improve performance by computing the links only in the specified interval
      - Add the integrated reaction matrix computation
      - Add the integrated graph for StoRtoS analysis 
      - TODO : same for StoS and RtoS analysis

    """

    logger.info("Performing analysis on the path")
    logger.info("-------------------------------")

    # Setting parameters if not user defined
    if t_init is None:
        t_init = [0, 0]
    if t_end is None:
        t_end = [1e4, 1e4]
    if spec2plot is None:
        spec2plot = ["all"]

    if analys == 'StoS':
        logger.info("Analysis -> Species to Species flux")
    elif analys == 'RtoS':
        logger.info("Analysis -> Reaction to Species flux")
    elif analys == 'StoRtoS':
        logger.info("Analysis -> Species to Reaction to Species flux")
    elif analys == 'Glob':
        logger.info("Analysis -> Global representation of the kinematic scheme")
    else:
        logger.error("ERROR : This analysis is not specified, please choose between StoS, RtoS or StoRtoS !")

    if method == 'PFA':
        logger.info('Method   -> Path Flux Analysis')
    elif method == 'ATOM':
        logger.info('Method   -> Atom Flux Analysis')
        if analys == "RtoS" or analys == "StoRtoS":
            logger.error("ERROR : The method ATOM cannot be used with this type of analys !")
            logger.error("==> Use the keyword : method = 'PFA'. ")
            exit()
        # Mettre le test de la méthode RtoS et StoRtoS
    else:
        logger.error("ERROR : This type of method is not available. "
                     "Please choose PFA (Path Flux Analysis) or ATOM (Atom Flux Analysis)")

    # Create the folder where the image + video are gonna be put in
    # Create the folder where the sensibility restore file is gonna be put in
    if not directory:
        directory = analysis_dir + '/analysis_graphs/'
    database.create_dir(directory, overwrite=False)

    # Initializations
    images = []
    sdbss = []
    maxim = 0
    old_grid = 0.0

    importance = 0.4
    # Check size of the plot
    if size == 'QHD':
        size_vector = [2560, 1440]
    elif size == 'FullHD':
        size_vector = [1920, 1080]
    elif size == 'HD':
        size_vector = [1280, 720]
    elif size == 'qHD':
        size_vector = [960, 540]
    elif size == 'YouTube144p':
        size_vector = [256, 144]
    elif isinstance(size, list):
        size_vector = [size[0], size[1]]
    else:
        logger.error('ERROR ! You should specify either the size of the vector or the keywords '
                     'QHD, FullHD, HD, qHD, Youtube144p.')
        sys.exit()

    # Create sample database on cases and set-up solution file
    logger.info("=> Creating sampling database")
    sdbs = samples_database(case_list, mechanism, integrated)
    beg = 0

    # Build the database of databases to compare 2 cases
    # Maybe add a loop to calculate several mechanisms
    for nb_case, case in enumerate(case_list):
        sdbss.append(sdbs[beg:beg + case.nsamples])
        beg += case.nsamples
    ct_mechanism = mechanism.ctmech

    # Loop on sdbss to analyse each case
    for index_sdb, sdb in enumerate(sdbss):

        # Time / Position used to compute the mean graph
        tstart = t_init[index_sdb]
        tend = t_end[index_sdb]

        # Create dictionary of species
        dic_label_spec = {}
        for indx in range(0, mechanism.ns):
            dic_label_spec[indx] = str(ct_mechanism.species_name(indx))

        # Create dictionary of reactions
        dic_label_reac = {}
        for indx in range(0, mechanism.nr):
            dic_label_reac[indx] = str(ct_mechanism.reaction(indx))

        # Get the maximum of the DRG fluxes
        if analys == 'StoS' or analys == 'Glob':
            maxim = get_max(analys, method, mechanism, sdb, maxim, dic_label_spec, spec, atom, tstart, tend)
        else:
            maxim = get_max(analys, method, mechanism, sdb, maxim, dic_label_reac, spec, tstart=tstart, tend=tend)

        # Compute linking matrix
        logger.info("=> Computing link matrix")
        first_int = True
        for i_sample, my_sample in enumerate(sdb[0:-1]):
            if tstart <= my_sample.grid <= tend:
                logger.info("=> Grid point is " + '{0:.10e}'.format(my_sample.grid))
                if method == 'PFA':
                    matrix_spec, matrix_reac = pfa_path(my_sample, mechanism)
                elif method == 'ATOM':
                    matrix_spec, matrix_reac = atom_path(my_sample, mechanism, atom)
                else:
                    matrix_spec = []
                    matrix_reac = []
                    logger.error('Your method is not implemented for the moment !')

                # Compute integral matrix
                # if i_sample == 0:
                if first_int:
                    # matrix_spec_int = matrix_spec * my_sample.grid
                    matrix_spec_int = matrix_spec * 0.0
                    # matrix_reac_int = matrix_reac * my_sample.grid
                    matrix_reac_int = matrix_reac * 0.0
                    old_grid = my_sample.grid
                    first_int = False
                else:
                    if tstart < my_sample.grid < tend:
                        matrix_spec_int += matrix_spec * (my_sample.grid - old_grid)
                        matrix_reac_int += matrix_reac * (my_sample.grid - old_grid)
                        old_grid = my_sample.grid
                    else:
                        old_grid = my_sample.grid

                # Construction of the graph

                # Saving graph for each sample

                # Create graph showing interactions
                digraph = Digraph(format=fmt)
                ratio = size_vector[1] / ((1 - importance) * size_vector[0])
                digraph.attr(ratio=str(ratio))
                logger.info("=> Drawing graph for links")
                if analys == 'StoRtoS':
                    digraph = reaction_production_graph(digraph, mechanism, matrix_reac, dic_label_reac, spec, maxim,
                                                        threshold, method)
                elif analys == 'RtoS':
                    digraph = reaction_graph(digraph, mechanism, matrix_reac, dic_label_reac, spec[0], maxim,
                                             threshold, method)
                elif analys == 'StoS':
                    digraph = production_graph(digraph, mechanism, matrix_spec, dic_label_spec, spec[0], maxim,
                                               threshold, method)
                elif analys == 'Glob':
                    digraph = global_graph(digraph, mechanism, matrix_spec, dic_label_spec, spec[0], maxim, threshold,
                                           method, atom, index_sdb)

                digraph.render(directory + '/graph', view=False, format=fmt)
                os.remove(directory + '/graph')

                logger.info("\n=> Drawing graph for case")
                # Create aside graph showing the physics
                ready_made_plot(directory, 'plot1.' + fmt, integrated, ct_mechanism, sdb, my_sample, importance,
                                size_vector, zoom=False)

                logger.info("=> Drawing graph for case (near view)")
                # Create aside zoom graph showing the physics
                ready_made_plot(directory, 'plot2.' + fmt, integrated, ct_mechanism,
                                sdb, my_sample, importance, size_vector, zoom=True)

                logger.info("=> Merging graphs")
                # Merge the two plot graphs (physic of the flame + zoom)
                x1, x2 = int(size_vector[0] * importance), int(size_vector[0] * importance)
                y1, y2 = int(size_vector[1] / 2), int(size_vector[1] / 2)
                img_merge = merge_vertical('plot1.' + fmt, x1, y1, 'plot2.' + fmt, x2, y2, directory)
                img_merge.save(directory + '/plot.' + fmt)

                # Create merge of the two graphs (previous merge + digraph)
                x1, x2 = int(size_vector[0] * (1 - importance)), int(size_vector[0] * importance)
                y1, y2 = int(size_vector[1]), int(size_vector[1])
                img_merge = merge_horizontal('graph.' + fmt, x1, y1, 'plot.' + fmt, x2, y2, directory)
                img_merge.save(directory + '/' + analys + '_' + '{:.10f}'.format(my_sample.grid) + '.' + fmt)
                images.append(analys + '_' + '{:.10f}'.format(my_sample.grid) + '.' + fmt)

        os.remove(directory + '/graph.' + fmt)
        os.remove(directory + '/plot.' + fmt)

        # Create the movie if asked
        if movie is True:
            from cv2 import VideoWriter, VideoWriter_fourcc, imread, resize

            logger.info("=> Creating movie")

            # Chose the folder to save the movie
            image_folder = directory
            video_name = directory + '/' + analys + '.avi'

            # Write the movie
            fourcc = VideoWriter_fourcc(*"XVID")
            video = VideoWriter(video_name, fourcc, float(1), (size_vector[0], size_vector[1]))
            for image in images:
                video.write(imread(os.path.join(image_folder, image)))

            # Save the movie
            video.release()

        # Create the integrated graph
        if analys == 'Glob':

            # Initialization of the graph
            if index_sdb == 0:
                Gtot = Digraph(format=fmt)
                ratio = size_vector[1] / ((1 - importance) * size_vector[0])
                Gtot.attr(ratio=str(ratio))
            logger.info("=> Drawing integrated graph for links")

            # Normalization with percentages
            for i in range(len(matrix_spec_int[:, 0])):
                sum_fluxes = np.sum(matrix_spec_int[i, :])
                if sum_fluxes > 0.0:
                    for j in range(len(matrix_spec_int[0, :])):
                        matrix_spec_int[i, j] = matrix_spec_int[i, j] / sum_fluxes * 100

            # Creation of the graph
            maxim = np.max(matrix_spec_int)
            Gtot = global_graph(Gtot, mechanism, matrix_spec_int, dic_label_spec, spec[0], 100, threshold, method, atom,
                                index_sdb, percentage=True, spec2plot=spec2plot)
            Gtot.render(directory + '/graph_integrated', view=False, format=fmt)

        if analys == 'StoRtoS':
            # Initialization of the graph
            if index_sdb == 0:
                Gtot = Digraph(format=fmt)
                ratio = size_vector[1] / ((1 - importance) * size_vector[0])
                Gtot.attr(ratio=str(ratio))
            logger.info("=> Drawing integrated graph for links")

            # # Normalization with percentages
            # for i in range(len(matrix_reac_int[:,0])):
            #     sum_fluxes = np.sum(matrix_reac_int[i,:])
            #     if sum_fluxes > 0.0 :
            #         for j in range(len(matrix_reac_int[0,:])):
            #             matrix_reac_int[i,j] = matrix_reac_int[i,j]/sum_fluxes*100

            # Creation of the graph
            maxim = np.max(matrix_reac_int)
            Gtot = reaction_production_graph(Gtot, mechanism, matrix_reac_int, dic_label_reac, spec, maxim, threshold,
                                             method)
            Gtot.render(directory + '/graph_integrated', view=False, format=fmt)

        logger.info("=> Done")


def pfa_path(my_sample, mechanism):
    """ Function to compute the coefficients of a Path Flux Analysis for
    path_analysis function.

    :param my_sample: sample object
    :param mechanism: mechanism object

    :return matrix_spec, matrix_reaction : path flux matrices for species (ns*ns+1) and reactions (nr*ns+1)
    (arrays of double)

    Created: 29/08/2019 [JW]

    """

    # Parameters
    my_ns = mechanism.ns
    my_nr = mechanism.nr
    net = mechanism.network
    ct_mechanism = mechanism.ctmech
    nup = net.nup
    nur = net.nur
    nu = nup - nur

    # Get reaction rates for each reaction
    ct_mechanism.TPY = float(my_sample.T), float(my_sample.P), [float(myY) for myY in my_sample.Y]
    omega = ct_mechanism.net_rates_of_progress

    # Get enthalpy from production and consumption of species in reactions
    hr_species = species_heat_release(mechanism, my_sample)
    i_hr = my_ns

    # Get production and consumption rates for species
    pa = ct_mechanism.creation_rates
    ca = ct_mechanism.destruction_rates

    # Workspace
    matrix_species = np.zeros([my_ns + 1, my_ns], 'd')
    matrix_reaction = np.zeros([my_ns, my_nr], 'd')

    # Evaluate modified DIC(i,j) of the species
    for i in range(my_ns):
        # reactions containing species i
        bool_i = net.indr[net.inds == i]
        ind_j = net.indsj[net.indsi == i]

        for j in ind_j:
            bool_j = net.indr[net.inds == j]  # reactions containing species j
            ind_k = np.intersect1d(bool_i, bool_j)  # reactions containing species i and j

            # Compute the DIC
            for k in ind_k:
                if np.sign(nu[i, k]) != np.sign(nu[j, k]):  # here is the modification
                    matrix_species[i, j] += (nup[i, k] - nur[i, k]) * omega[k]

            # Normalize
            matrix_species[i, j] = matrix_species[i, j] / max(max(pa[i], ca[i]), 1e-60)

    for i in range(my_ns):
        matrix_species[i_hr, i] += abs(hr_species[i])

    # Evaluate DIC(i,j) of the reactions
    for i in range(my_ns):
        # reactions containing species i
        ind_k = net.indr[net.inds == i]  # reactions containing species i

        # Compute the DIC
        for k in ind_k:
            matrix_reaction[i, k] = (nup[i, k] - nur[i, k]) * omega[k]

            # Normalize
            matrix_reaction[i, k] = matrix_reaction[i, k] / max(max(pa[i], ca[i]), 1e-60)

    return matrix_species, matrix_reaction


def atom_path(mysample, mechanism, atom='C'):
    """ Function to compute the coefficients of an Atom Flux Analysis for
    path_analysis function.

    :param mysample: sample object
    :param mechanism: mechanism object
    :param atom: allow to choose the atom in the Atom flux method

    :return DIC_spec, DIC_reac : atom flux matrices for species (ns*ns+1) and reactions (nr*ns+1) (arrays of double)

    Created: 29/08/2019 [JW]
    Last modification 06/01/2020 [AP]
    - Correction of ATOM flux computation
    - Typo correction to do (mysample, myns, mechanism, ctmech,...) Old syntax keeped to avoid issues
    - Comments
    """

    # Parameters
    myns = mechanism.ns
    mynr = mechanism.nr
    net = mechanism.network
    ctmech = mechanism.ctmech
    nup = net.nup
    nur = net.nur

    # atom = 'C'
    # Get reaction rates for each reaction
    ctmech.TPY = float(mysample.T), float(mysample.P), [float(myY) for myY in mysample.Y]
    omega = ctmech.net_rates_of_progress

    # Workspace
    DIC_spec = np.zeros([myns, myns], 'd')
    DIC_reac = np.zeros([myns + 1, mynr], 'd')

    # Evaluate DIC(i,j) of the species
    for i in range(myns):
        # reactions containing species i
        booli = net.indr[net.inds == i]
        indj = net.indsj[net.indsi == i]

        for j in indj:
            boolj = net.indr[net.inds == j]  # reactions containing species j
            indk = np.intersect1d(booli, boolj)  # reactions containing species i and j

            # Compute the reaction path A_ijk = (nA_i*nA_j*r_k)/NA_k
            for k in indk:
                N = 0
                for r in ctmech.reaction(k).reactants:
                    index = ctmech.species_index(r)
                    N += ctmech.n_atoms(r, atom) * nur[index, k]
                # Forward reactions
                if nur[i, k] != 0 and nup[j, k] != 0 and omega[k] > 0:
                    if N != 0:
                        DIC_spec[i, j] += ctmech.n_atoms(ctmech.species_name(i), atom) * nur[i, k] * \
                                          ctmech.n_atoms(ctmech.species_name(j), atom) * nup[j, k] * omega[k] / N
                # Backward reactions
                if nur[j, k] != 0 and nup[i, k] != 0 and omega[k] < 0:
                    if N != 0:
                        DIC_spec[i, j] -= ctmech.n_atoms(ctmech.species_name(i), atom) * nur[j, k] * \
                                          ctmech.n_atoms(ctmech.species_name(j), atom) * nup[i, k] * omega[k] / N

    return DIC_spec, DIC_reac


def get_max(analys, method, mechanism, sdb, maxim, dic_label, spec, atom='C', tstart=0.0, tend=10000.0):
    """ Function to compute the maximum coefficients of a given method.

    :param analys: analysis used (string)
    :param method: method used (string)
    :param mechanism: mechanism object
    :param sdb: sample database used
    :param maxim: maximum object (double)
    :param dic_label: dictionary label for naming reactions (string array)
    :param spec: species to consider (string)
    :param atom: atom in the Atom flux method
    :param tstart: Initial time or position for integrated graph (floats)
    :param tend: Final time or position for integrated graph (floats)

    :return maxim: maximum object (double)

    Created: 29/08/2019 [JW]
    Last modification 29/04/2020 [AP]
    - Add new parameter tstart and tend to get maximum only in the interval of interest

    """

    # Parameters
    ct_mechanism = mechanism.ctmech
    matrix_reaction = []
    matrix_species = []

    if analys == 'StoRtoS':
        for i_sample, my_sample in enumerate(sdb):
            if tstart <= my_sample.grid <= tend:
                if method == 'PFA':
                    matrix_species, matrix_reaction = pfa_path(my_sample, mechanism)
                elif method == 'ATOM':
                    matrix_species, matrix_reaction = atom_path(my_sample, mechanism)
                for s0 in spec:
                    ind_e = ct_mechanism.species_index(s0)
                    for s1 in spec:
                        for ind_f, f in enumerate(dic_label):
                            if ((s1 in ct_mechanism.reaction(ind_f).reactants and
                                 s0 in ct_mechanism.reaction(ind_f).products)
                                    or (s0 in ct_mechanism.reaction(ind_f).reactants and
                                        s1 in ct_mechanism.reaction(ind_f).products)):
                                maxim = np.maximum(maxim, abs(matrix_reaction[ind_e, ind_f]))

    elif analys == 'RtoS':
        for i_sample, my_sample in enumerate(sdb):
            if tstart <= my_sample.grid <= tend:
                if method == 'PFA':
                    matrix_species, matrix_reaction = pfa_path(my_sample, mechanism)
                elif method == 'ATOM':
                    matrix_species, matrix_reaction = atom_path(my_sample, mechanism)
                for s0 in spec:
                    ind_e = ct_mechanism.species_index(s0)
                    for ind_f, f in enumerate(dic_label):
                        if matrix_reaction[ind_e, ind_f] != 0:
                            maxim = np.maximum(maxim, abs(matrix_reaction[ind_e, ind_f]))

    elif analys == 'StoS' or analys == 'Glob':
        for i_sample, my_sample in enumerate(sdb):
            if tstart <= my_sample.grid <= tend:
                if method == 'PFA':
                    matrix_species, matrix_reaction = pfa_path(my_sample, mechanism)
                elif method == 'ATOM':
                    matrix_species, matrix_reaction = atom_path(my_sample, mechanism, atom)
                for s0 in spec:
                    ind_e = ct_mechanism.species_index(s0)
                    for ind_f, f in enumerate(dic_label):
                        if matrix_species[ind_e, ind_f] != 0:
                            maxim = np.maximum(maxim, abs(matrix_species[ind_e, ind_f]))

    return maxim


def production_graph(digraph, mechanism, matrix_species, dic_label, spec, maxim, threshold, method='PFA'):
    """ Function to compute graph for StoS method.

    :param digraph: graph object (Digraph)
    :param mechanism: mechanism object
    :param matrix_species: matrix with links (array double)
    :param dic_label: dictionary label for naming reactions (string array)
    :param spec: species to consider (string)
    :param maxim: maximum object (double)
    :param threshold: threshold applied on fluxes
    :param method: method used, either PFA (path flux analysis) or ATOM (atom flux analysis) (string)

    :return digraph: graph object (Digraph)

    Created: 29/08/2019 [JW]
    Last modification 19/12/2019 [AP]
    - New parameter "method" to choose the Atom flux method
    """

    ct_mechanism = mechanism.ctmech
    ind_e = ct_mechanism.species_index(spec)
    e = ind_e
    digraph.attr('node', shape='doublecircle')
    eq_e = str(ct_mechanism.species_name(e))
    digraph.node(eq_e)

    # Allows to keep the right way when the arrows are displayed
    if method == 'ATOM':
        matrix_species = -matrix_species

    for ind_f, f in enumerate(dic_label):
        if matrix_species[ind_e, ind_f] != 0 and abs(matrix_species[ind_e, ind_f] / maxim) > threshold:

            if ct_mechanism.species_name(f) != ct_mechanism.species_name(e):
                digraph.attr('node', shape='circle')
            else:
                digraph.attr('node', shape='doublecircle')
            eq_f = str(ct_mechanism.species_name(f))
            digraph.node(eq_f)
            digraph.attr('node', shape='box')
            eq_w = 'w = ' + '{0:.2e}'.format(abs(matrix_species[ind_e, ind_f]))

            if abs(matrix_species[ind_e, ind_f]) / maxim > 0.1:
                color = "0.000 1.000 1.000 " + str(abs(matrix_species[ind_e, ind_f]) / maxim)
            else:
                color = "0.000 1.000 0.500 " + str(abs(matrix_species[ind_e, ind_f]) / maxim)

            if matrix_species[ind_e, ind_f] < 0:
                digraph.edge(eq_e, eq_w, penwidth=str(-1 / (np.log10(-matrix_species[ind_e, ind_f] / maxim) - 1) * 3),
                             color=color)
                digraph.edge(eq_w, eq_f, penwidth=str(-1 / (np.log10(-matrix_species[ind_e, ind_f] / maxim) - 1) * 3),
                             color=color)
            else:
                digraph.edge(eq_f, eq_w, penwidth=str(-1 / (np.log10(matrix_species[ind_e, ind_f] / maxim) - 1) * 3),
                             color=color)
                digraph.edge(eq_w, eq_e, penwidth=str(-1 / (np.log10(matrix_species[ind_e, ind_f] / maxim) - 1) * 3),
                             color=color)

    return digraph


def reaction_graph(digraph, mechanism, matrix_reaction, dic_label, spec, maxim, threshold, method='PFA'):
    """ Function to compute graph for StoS method.

    :param digraph: graph object (Digraph)
    :param mechanism: mechanism object
    :param matrix_reaction: matrix with links (array double)
    :param dic_label: dictionary label for naming reactions (string array)
    :param spec: species to consider (string)
    :param maxim: maximum object (double)
    :param threshold: threshold applied on fluxes
    :param method: method used, either PFA (path flux analysis) or ATOM (atom flux analysis) (string) ! TO REMOVE !)

    :return digraph: graph object (Digraph)

    Created: 29/08/2019 [JW]
    Last modification 19/12/2019 [AP]
    - New parameter "method" to choose the Atom flux method
    """

    ct_mechanism = mechanism.ctmech
    ind_e = ct_mechanism.species_index(spec)
    e = ind_e
    digraph.attr('node', shape='doublecircle')
    digraph.node(str(ct_mechanism.species_name(e)))

    # Allows to keep the right way when the arrows are displayed
    # /!\ TO REMOVE atom flux method can't be used with reactions /!\
    if method == 'ATOM':
        matrix_reaction = -matrix_reaction

    for ind_f, f in enumerate(dic_label):
        if matrix_reaction[ind_e, ind_f] != 0 and abs(matrix_reaction[ind_e, ind_f] / maxim) > threshold:

            digraph.attr('node', shape='box')
            eq = str(ct_mechanism.reaction(f)) + '\n w = ' + '{0:.2e}'.format(abs(matrix_reaction[ind_e, ind_f]))
            digraph.node(eq)

            if abs(matrix_reaction[ind_e, ind_f]) / maxim > 0.1:
                color = "0.000 1.000 1.000 " + str(abs(matrix_reaction[ind_e, ind_f]) / maxim)
            else:
                color = "0.000 1.000 0.500 " + str(abs(matrix_reaction[ind_e, ind_f]) / maxim)

            if matrix_reaction[ind_e, ind_f] < 0:
                digraph.edge(str(ct_mechanism.species_name(e)), str(eq),
                             penwidth=str(-1 / (np.log10(-matrix_reaction[ind_e, ind_f] / maxim) - 1)), color=color)
            else:
                digraph.edge(eq, str(ct_mechanism.species_name(e)),
                             penwidth=str(-1 / (np.log10(matrix_reaction[ind_e, ind_f] / maxim) - 1)), color=color)

    return digraph


def reaction_production_graph(digraph, mechanism, matrix_reaction, dic_label, spec, maxim, threshold, method='PFA'):
    """ Function to compute graph for StoS method.

    :param digraph: graph object (Digraph)
    :param mechanism: mechanism object
    :param matrix_reaction: matrix with links (array double)
    :param dic_label: dictionary label for naming reactions (string array)
    :param spec: species to consider (string)
    :param maxim: maximum object (double)
    :param threshold: threshold applied on fluxes
    :param method: method used, either PFA (path flux analysis) or ATOM (atom flux analysis) (string) ! TO REMOVE !

    :return digraph: graph object (Digraph)

    Created: 29/08/2019 [JW]
    Last modification 29/04/2020 [AP]
    - New input parameters to be able to chose the maximum and the threashold
    """

    # maxim = 1
    # threshold = 0
    net = mechanism.network
    ct_mechanism = mechanism.ctmech

    # Allows to keep the right way when the arrows are displayed
    # /!\ TO REMOVE atom flux method can't be used with reactions /!\
    if method == 'ATOM':
        matrix_reaction = -matrix_reaction

    for s0 in spec:
        ind_e = ct_mechanism.species_index(s0)
        e = ind_e
        digraph.attr('node', shape='doublecircle')
        digraph.node(str(ct_mechanism.species_name(e)))

        for s1 in spec:
            ind_e_2 = ct_mechanism.species_index(s1)
            if s0 is not s1:
                for ind_f, f in enumerate(dic_label):
                    if net.deltaSR[ind_e][ind_f] != 0 and net.deltaSR[ind_e_2][ind_f] != 0 \
                            and matrix_reaction[ind_e, ind_f] != 0 \
                            and abs(matrix_reaction[ind_e, ind_f] / maxim) > threshold:
                        digraph.attr('node', shape='box')
                        eq = str(ct_mechanism.reaction(f))
                        digraph.node(eq)

                        if abs(matrix_reaction[ind_e, ind_f]) / maxim > 0.1:
                            color = "0.000 1.000 1.000 " + str(abs(matrix_reaction[ind_e, ind_f]) / maxim)
                        else:
                            color = "0.000 1.000 0.500 " + str(abs(matrix_reaction[ind_e, ind_f]) / maxim)

                        if matrix_reaction[ind_e, ind_f] < 0:
                            digraph.edge(str(ct_mechanism.species_name(e)), str(eq),
                                         penwidth=str(-1 / (np.log10(-matrix_reaction[ind_e, ind_f] / maxim) - 1)),
                                         color=color,
                                         label='w = ' + '{0:.2e}'.format(abs(matrix_reaction[ind_e, ind_f])))
                        else:
                            digraph.edge(eq, str(ct_mechanism.species_name(e)),
                                         penwidth=str(-1 / (np.log10(matrix_reaction[ind_e, ind_f] / maxim) - 1)),
                                         color=color,
                                         label='w = ' + '{0:.2e}'.format(abs(-matrix_reaction[ind_e, ind_f])))

    return digraph


def global_graph(G, mechanism, DIC_spec, dic_label, spec='CH4', maxim=1, seuil=1e-1, method='ATOM', atom='C',
                 index_sdb=0,
                 percentage=False, spec2plot=None):
    """ Function to build the graph of the kinetic scheme.

    :param G: graph object (Digraph)
    :param mechanism: mechanism object
    :param DIC_spec: matrix with links (array double)
    :param dic_label: dictionary label for naming reactions (string array)
    :param spec: species to consider (string)
    :param maxim: maximum object (double)
    :param seuil: threshold applied on fluxes
    :param method: method used, either PFA (path flux analysis) or ATOM (atom flux analysis) (string)
    :param atom: atom used for the ATOM method (string)
    :param index_sdb: index used for the color of the arrows (boolean)
    :param percentage: activate or not a percentage analysis
    :param spec2plot: Species to display in the graph (string array)

    :return G: graph object (Digraph)

    Created: 01/10/2019 [AP]
    Last modification 29/04/2020 [AP]
    - New input parameter spec2plot to choose the species to display in the global graph
    The value is set to 'all' when the graph correspond to a unique sample, but the species can be
    chosen for the integrated graph.
    - Add new logger information on the sum of the path coefficients to be able to know if the result 
    displayed if more of less close to the reality. (in case of percentage graph, the sum of the coefficient should be
    equal to 100)

    """

    # Extract mechanism
    ctmech = mechanism.ctmech

    # Initialize spec2plot
    if spec2plot is None:
        spec2plot = ["all"]

    # Legend
    G.attr('node', shape='box')
    if percentage:
        msg1 = str("Atom flux of " + atom + "\n Integrated")
    else:
        msg1 = str("Atom flux of " + atom + "\n T = " + str(int(ctmech.T)) + ' K')
        spec2plot = "all"
    G.node('L', msg1)

    # Create the first node of the graph with the species of interest
    ind_e = ctmech.species_index(spec)
    e = ind_e
    G.attr('node', shape='doublecircle')
    eq_e = str(ctmech.species_name(e))
    G.node(eq_e)
    G.edge('L', eq_e, penwidth="0", arrowhead="none", color="white")

    # Allows to keep the right way when the arrows are displayed
    if method == 'ATOM':
        DIC_spec = -DIC_spec

    # sub_spec is the list with the new species at each step
    sub_spec = []
    # tree is the list with the species which we want to study the atom flux
    tree = [spec]
    # tot_spec is the list with all the species already studied
    tot_spec = tree

    # Beginning of the loop to build the tree
    stop = False
    generation = 0
    while not stop:
        logger.info("\n# Generation : " + str(generation))

        # loop on the species which we want to study the atom flux
        for eq_e in tree:
            if "all" in spec2plot or eq_e in spec2plot:

                ind_e = ctmech.species_index(eq_e)

                # loop on all the species to determine all the influences
                sum_displayed = 0
                for ind_f, f in enumerate(dic_label):
                    if "all" in spec2plot or str(ctmech.species_name(f)) in spec2plot:

                        # A link is created in the tree only above a specific threashold
                        if DIC_spec[ind_e, ind_f] != 0 and abs(DIC_spec[ind_e, ind_f] / maxim) > seuil:

                            # Define the new species to study at the next "generation"
                            if ctmech.species_name(f) not in tot_spec and ctmech.species_name(f) not in sub_spec:
                                sub_spec.append(ctmech.species_name(f))

                            # Create a node for new species
                            if ctmech.species_name(f) != ctmech.species_name(e):
                                G.attr('node', shape='circle')
                            else:
                                G.attr('node', shape='doublecircle')
                            eq_f = str(ctmech.species_name(f))
                            if eq_f not in tot_spec:
                                G.node(eq_f)

                            # Set the color of the arrow between species
                            if index_sdb == 0:
                                color = "0.000 1.000 1.000 " + str(abs(DIC_spec[ind_e, ind_f]) / maxim)
                            else:
                                color = "0.000 0.000 0.000 " + str(abs(DIC_spec[ind_e, ind_f]) / maxim)

                            # print information in the logger
                            logger.info("Link " + str(eq_e) + "==>" + str(eq_f) + " : ")
                            logger.info("coeff :" + str(-DIC_spec[ind_e, ind_f]))
                            sum_displayed += -DIC_spec[ind_e, ind_f]

                            # Create the link between the species
                            if str(eq_e) not in ["CO2", "H2O"]:  # Otherwise it can lead to infinite loops
                                # O2 ==> CO ==> CO2 ==> CO ...
                                # if str(eq_e)!="CCC":
                                if DIC_spec[ind_e, ind_f] < 0:
                                    if percentage:
                                        G.edge(eq_e, eq_f, label='{:.1f}'.format(abs(DIC_spec[ind_e, ind_f])) + " %",
                                               penwidth="1", color=color)
                                    else:
                                        G.edge(eq_e, eq_f, label='{:.1f}'.format(abs(DIC_spec[ind_e, ind_f])),
                                               penwidth=str(abs(DIC_spec[ind_e, ind_f]) / maxim * 5), color=color)
                                else:
                                    if percentage:
                                        G.edge(eq_f, eq_e, label='{:.1f}'.format(abs(DIC_spec[ind_e, ind_f])) + " %",
                                               penwidth="1", color=color)
                                    else:
                                        G.edge(eq_f, eq_e, label='{:.1f}'.format(abs(DIC_spec[ind_e, ind_f])),
                                               penwidth=str(abs(DIC_spec[ind_e, ind_f]) / maxim * 5), color=color)
                if percentage:
                    logger.info("Missing fluxes for " + str(eq_e) + " : " + str(100 - sum_displayed) + " %")
                else:
                    logger.info("Sum fluxes for " + str(eq_e) + " : " + str(sum_displayed))

        # If no new species are found the loop is stopped
        if not sub_spec:
            stop = True

        # Otherwise a new set of species is defined for the next generation
        else:
            generation += 1
            for el in sub_spec:
                tot_spec.append(el)
            tree = sub_spec
            sub_spec = []

        # Safety condition to avoid infinite loop
        if generation > 10:
            stop = True

    return G


def merge_horizontal(image1, x1, y1, image2, x2, y2, directory):
    """ Function to merge a picture horizontally.

    :param image1: first image object
    :param x1: first dimension of first picture
    :param y1: second dimension of first picture
    :param image2: second image object
    :param x2: first dimension of second picture
    :param y2: second dimension of second picture

    :param directory: directory to store the picture in (string)

    :return img_merge: image object

    Created: 29/08/2019 [JW]

    """

    images_list = [directory + '/' + image1, directory + '/' + image2]
    img_open = [Image.open(i) for i in images_list]
    img_merge = np.hstack((np.asarray(img_open[0].resize([x1, y1], Image.ANTIALIAS)),
                           np.asarray(img_open[1].resize([x2, y2], Image.ANTIALIAS))))
    img_merge = Image.fromarray(img_merge)

    return img_merge


def merge_vertical(image1, x1, y1, image2, x2, y2, directory):
    """ Function to merge a picture vertically.

    :param image1: first image object
    :param x1: first dimension of first picture
    :param y1: second dimension of first picture
    :param image2: second image object
    :param x2: first dimension of second picture
    :param y2: second dimension of second picture

    :param directory: directory to store the picture in (string)

    :return img_merge: image object

    Created: 29/08/2019 [JW]

    """

    images_list = [directory + '/' + image1, directory + '/' + image2]
    img_open = [Image.open(i) for i in images_list]
    img_merge = np.vstack((np.asarray(img_open[0].resize([x1, y1], Image.ANTIALIAS)),
                           np.asarray(img_open[1].resize([x2, y2], Image.ANTIALIAS))))
    img_merge = Image.fromarray(img_merge)

    return img_merge


def ready_made_plot(directory, image, integrated, ct_mechanism, sdb, my_sample, importance, size_vector, zoom=False):
    """ Function to compute a ready-made plot.

    :param directory: directory to store the picture in (string)
    :param image: image used to plot the figure
    :param integrated: if the sample database is integrated or not (boolean)
    :param ct_mechanism: object solution of the mechanism
    :param sdb: sample database object
    :param my_sample: sample considered
    :param importance: float representing the importance of a graph compared to
    another in the vertical direction
    :param size_vector: float representing the sizes for the different resolutions
    :param zoom: boolean that determines if the field should be zoomed or not

    :return None

    Created: 29/08/2019 [JW]
    Last modification 25/06/2020 [AP]
    - Set the ylim to [-1;1] in case of strong endothermic step

    """

    # Get interesting values for the plot to get the dynamic of the flame
    grid = [m.grid for m in sdb]
    temperature = [m.T for m in sdb]
    mass_fraction_oxygen = [m.Y[ct_mechanism.species_index('O2')] for m in sdb]
    hr = [m.HR for m in sdb]

    # Get the values of the x-axis to zoom in the reaction zone
    i = 0
    if hr is not None:
        f = (hr - np.min(hr)) / (np.max(hr) - np.min(hr))
        df = np.gradient(f, grid)
        while df[i] < 1000 and i < len(df) - 1:
            i = i + 1
        if i < len(df) - 1:
            x_zoom_1 = grid[i]
        else:
            x_zoom_1 = min(grid)
        i = len(sdb) - 1
        while df[i] > -1000 and i > 1:
            i = i - 1
        if i > 1:
            x_zoom_2 = grid[i]
        else:
            x_zoom_2 = max(grid)
    else:
        x_zoom_1 = min(grid)
        x_zoom_2 = max(grid)

    # Plot the graph and save it
    plt.figure(figsize=[importance * 20, size_vector[1] / size_vector[0] * 10])
    plt.plot(grid, temperature / np.max(temperature), '+-')
    plt.plot(grid, mass_fraction_oxygen / np.max(mass_fraction_oxygen), '+-')
    plt.plot(grid, hr / np.max(hr), '+-')
    plt.xlabel('Grid (t or x)', fontsize=15)
    plt.ylabel('Dimensionless data', fontsize=15)
    plt.xlim([x_zoom_1, x_zoom_2])
    plt.legend(['Temperature', 'Oxygen Mass Fraction', 'Heat Release'], fontsize=15)
    # plt.ylim([-1, 1])
    if zoom is True:
        plt.xlim([x_zoom_1, x_zoom_2])
        plt.ylim([-1, 1])
        if x_zoom_1 < my_sample.grid < x_zoom_2:
            plt.axvline(x=my_sample.grid)
    else:
        plt.xlim([np.min(grid), np.max(grid)])
        plt.ylim([-1, 1])
        plt.axvline(x=my_sample.grid)
    if integrated is True:
        plt.title('Integrated characteristics of the flame', fontsize=20)
        if zoom is True:
            plt.title('Integrated characteristics of the flame (zoom)', fontsize=20)
    else:
        plt.title('Characteristics of the flame : x = ' + '{:.2e}'.format(my_sample.grid), fontsize=20)
        if zoom is True:
            plt.title('Characteristics of the flame (zoom) : x = ' + '{:.2e}'.format(my_sample.grid), fontsize=20)

    plt.savefig(directory + "/" + image)
    plt.close()


###### Rests of the old sampling process to be replaced in this module #####

class Sample:
    """
    All methods and attributes of a chemical sample, as used in chemistry reduction
    - TODO : remove this class so that analysis is working with the sampling.py file
    """

    def __init__(self, case, mechanism, Y, P, T, HR, grid=None):
        """
        Generator of "sample" objects. A sample is defined as a container with:
        * a cantera Solution object according to which compositions are obtained
        * a composition state, with a mass fractions vector, a temperature and a pressure

        :param case: Case object
        :param mechanism: Mechanism object
        :param Y: mass fractions
        :param P: pressure
        :param T: temperature
        :param HR: heat release rate
        :param grid: grid
        :param grid: id of the case the sample belongs to

        Created: 11/14/17 [PP]
        Last modified: 11/14/17 [PP]
        """

        self.case = case
        self.mechanism = mechanism
        self.Y = Y
        self.T = T
        self.P = P
        self.HR = HR
        self.grid = grid


def samples_database(caselist, mechanism, integrated=True):
    """
    Go through all case instances and extract chemical samples from each of them.
    Results are concatenated into common database

    :param caselist: list of Case objects
    :param mechanism: Mechanism object containing all info about mechanism
    :param integrated: if False, takes instantaneous data

    :return sdb: combined sample database

    Created: 17/11/14 [PP]
    Last modified: 17/11/14 [PP]
    """

    sdb = []
    for case in caselist:
        sdb += extract_samples(case, mechanism, integrated)

    return sdb


def scalars_database(caselist, mechanism):
    """
    Go through all case instances and extract chemical scalars from each of them.
    Results are concatenated into common database

    :param caselist: list of Case objects
    :param mechanism: Mechanism object containing all info about mechanism

    :return sdb: combined sample database

    Created: 17/11/14 [PP]
    Last modified: 17/11/14 [PP]
    """

    scalars_db = np.zeros([len(caselist), mechanism.ctmech.n_species])
    for index_case, case in enumerate(caselist):
        data_dict = case.data_dict(mechanism)
        for index_spec, spec in enumerate(mechanism.species_names):
            scalars_db[index_case, index_spec] = tools.extract_scalar(data_dict[spec], 'int', grid=data_dict['grid'])

    return scalars_db


def extract_samples(case, mechanism, integrated=True):
    """
    Extracts samples for dumped profiles

    :param case: Case object
    :param mechanism: Mechanism object containing all info about mechanism
    :param integrated: if False, takes instantaneous data

    :return: collection of samples
    """
    # Initializing samples collection array
    samples_collection = []

    if not hasattr(case, 'ndata'):
        case.run(mechanism)

    data = case.extract_profile_old(mechanism)
    data_names_dict = case.names_dictionary_old(mechanism)

    # Checking if the data is consistent
    all_species_in_data = True
    for spec in mechanism.species_names:
        if spec not in data_names_dict:
            all_species_in_data = False

    if not all_species_in_data:
        case.run(mechanism, overwrite=True)
        data = case.extract_profile_old(mechanism)
        data_names_dict = case.names_dictionary_old(mechanism)

    # Skipping first result line as the initialization values can cause errors
    data = data[1:, :]
    case.ndata -= 1

    n = 0
    ncur = 0
    nstart = ncur
    nspecies = mechanism.ns
    axis = 0

    q, r = divmod(case.ndata, case.nsamples)

    if integrated:

        grid_array = []
        T_array = []
        P_array = []
        Y_array = []
        HR_array = []

        # Loop through the rest
        for line in data:
            # Clip very small concentrations
            line[1:] = [0.0 if abs(float(item)) < 1e-15 else item for item in line[1:]]

            # Arrays for integration
            grid_array.append(float(line[data_names_dict['Grid']]))
            T_array.append(float(line[data_names_dict['Temperature']]))
            P_array.append(float(line[data_names_dict['Pressure']]))
            Y_array.append([float(line[data_names_dict[spec]]) for spec in mechanism.species_names])
            if 'HeatRelease' in data_names_dict:
                HR_array.append(float(line[data_names_dict['HeatRelease']]))
            else:
                HR_array = []

            # Handle fact that number of data is not divisible by q
            if ncur - nstart < r:
                myq = q + 1
            else:
                myq = q

            if n % myq == 0 and n not in [0, 1]:
                if len(grid_array) == 1:
                    P = P_array[0]
                    T = T_array[0]
                    Y = Y_array[0]
                    if 'HeatRelease' in data_names_dict:
                        HR = HR_array[0]
                    else:
                        HR = 0
                else:
                    interval = grid_array[-1] - grid_array[0]
                    axis = (grid_array[-1] + grid_array[0]) / 2
                    P = np.trapz(P_array, grid_array) / interval
                    T = np.trapz(T_array, grid_array) / interval
                    if 'HeatRelease' in data_names_dict:
                        HR = np.trapz(HR_array, grid_array) / interval
                    else:
                        HR = 0
                    Y_array = np.transpose(Y_array)
                    Y = [np.trapz(Yspec_array, grid_array) / interval for Yspec_array in Y_array]

                # Reset the arrays
                grid_array = []
                P_array = []
                T_array = []
                HR_array = []
                Y_array = []

                mysample = Sample(case=case, mechanism=mechanism,
                                  Y=Y, P=P, T=T, HR=HR, grid=axis)

                samples_collection.extend([mysample])
                ncur += 1
            n += 1
    else:
        for line in data:
            # Clip very small concentrations
            line[1:] = [0.0 if abs(float(item)) < 1e-15 else item for item in line[1:]]
            # Handle fact that number of data is not divisible by q
            if ncur - nstart < r:
                myq = q + 1
            else:
                myq = q

            if n % myq == 0 and n not in [0, 1]:
                if 'HeatRelease' in data_names_dict:
                    HR = line[nspecies + 3]
                else:
                    HR = []
                mysample = Sample(case=case, mechanism=mechanism,
                                  Y=line[3:nspecies + 3], P=line[2], T=line[1],
                                  HR=HR, grid=line[0])

                samples_collection.extend([mysample])
                ncur += 1
            n += 1

    return samples_collection


def species_heat_release(mechanism, my_sample):
    """
    Computes the heat release rate from each species

    :param mechanism: Mechanism object
    :param sample: Sample object

    :return hr_species: list of heat release rate of each species
    """

    ctmech = mechanism.ctmech
    net = mechanism.network

    species_names = ctmech.species_names

    ctmech.TPY = float(my_sample.T), float(my_sample.P), [float(myY) for myY in my_sample.Y]

    hr_species = np.zeros((len(species_names)))
    hr_reactions = reactions_heat_release(mechanism, my_sample)
    for i_spec, spec in enumerate(species_names):
        # reactions containing species i
        booli = net.indr[net.inds == i_spec]
        for i_reac in booli:
            hr_species[i_spec] = hr_species[i_spec] + hr_reactions[i_reac]

    return hr_species


def reactions_heat_release(mechanism, my_sample):
    """
    Computes the heat release rate from each reaction

    :param mechanism: Mechanism object
    :param sample: Sample object

    :return hr_reactions: list of heat release rate of each reaction
    """

    ctmech = mechanism.ctmech

    n_reac = len(ctmech.reactions())

    ctmech.TPY = float(my_sample.T), float(my_sample.P), [float(myY) for myY in my_sample.Y]

    hr_reactions = [- ctmech.net_rates_of_progress[i_reac] * ctmech.delta_enthalpy[i_reac] for i_reac in range(n_reac)]

    return hr_reactions
