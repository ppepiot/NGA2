"""Module containing all qss-related functions"""

import ARCANE.display as display
import ARCANE.tools as tools
import networkx as nx
import numpy as np

logger = display.Logger()
logger.set_log('logCompute')


def print_constants_qss(print_variables, constants, qss_variables):
    """Prints the constants declaration specific to qss handling in the f90 file

    Parameters
    ----------
    print_variables : dict
        dictionary containing the f90 file to write in and code specific infos
    constants : dict
        dictionary of constants needed for computation
    qss_variables : dict
        dictionary of QSS related parameters needed for computation

    """

    f = print_variables['f']
    precision = print_variables['precision']

    nqss = constants['nqss']

    species_qss_names = qss_variables['species_qss_names']
    index_coupled_qss_i = qss_variables['index_coupled_qss_i']
    index_coupled_qss_j = qss_variables['index_coupled_qss_j']

    coupled_qss_matrix = qss_variables['coupled_qss_matrix']

    # Evaluation of QSS concentrations

    text = """\
  ! ----------------------------------------------- !
  ! Evaluation of QSS concentrations                !
  ! ----------------------------------------------- !
  subroutine get_qss(cqss,c,k,M)
    implicit none

    real({0}), dimension(nqss) :: cqss
    real({0}), dimension(nspec) :: c
    real({0}), dimension(nreac + nreac_reverse) :: k
    real({0}), dimension(nTB + nFO) :: M

    integer :: index

"""

    for i in range(nqss):
        text += '    real(' + precision + ') :: ' + species_qss_names[i] + '_ct' '\n'
        text += '    real(' + precision + ') :: ' + species_qss_names[i] + '_num' '\n'
        text += '    real(' + precision + ') :: ' + species_qss_names[i] + '_denom' '\n'

        for j in index_coupled_qss_j[index_coupled_qss_i == i]:
            if species_qss_names[i] != species_qss_names[j] and coupled_qss_matrix[i, j] == 1:
                text += '    real(' + precision + ') :: ' + species_qss_names[i] + '_' + species_qss_names[j] + '\n'
    text += '  \n'

    text = text.format(precision)

    f.write(text)


def qss_coupling(constants, mech_variables, reactions_variables, qss_variables, force=True):
    """Computes the need dependencies between QSS species

    Parameters
    ----------
    constants : dict
        dictionary constants needed for computation
    reactions_variables : dict
        dictionary of reactions related parameters needed for computation
    mech_variables : dict
        dictionary of mechanism related parameters needed for computation
    qss_variables : dict
        dictionary of QSS related parameters needed for computation
    force : bool
        if True, species and reactions will be removed from the mechanism if they induce QSS error (Default value = True)

    Returns
    -------
    coupled_qss_matrix
        a matrix with the value 1 if i needs j

    """

    nqss = constants['nqss']

    nur = mech_variables['nur']
    nup = mech_variables['nup']

    is_reversible = reactions_variables['is_reversible']

    species_qss_names = qss_variables['species_qss_names']
    qss_index_global = qss_variables['qss_index_global']
    index_coupled_qss_i = qss_variables['index_coupled_qss_i']
    index_coupled_qss_j = qss_variables['index_coupled_qss_j']
    index_qss_reactions_i = qss_variables['index_qss_reactions_i']
    index_qss_reactions_j = qss_variables['index_qss_reactions_j']

    coupled_qss_matrix = np.zeros([nqss, nqss])

    for i in range(nqss):
        token = []
        for j in index_coupled_qss_j[index_coupled_qss_i == i]:
            if j != i:
                for k in index_qss_reactions_j[index_qss_reactions_i == i]:
                    if ((nup[qss_index_global[i], k] >= 1 and nur[qss_index_global[j], k] >= 1) or
                        ((nur[qss_index_global[i], k] >= 1 and nup[qss_index_global[j], k] >= 1) and
                         is_reversible[k] == 1)) and (species_qss_names[j] not in token):
                        coupled_qss_matrix[i, j] = 1

    # Checks if quadratic coupling is involved
    exit_error = not force
    check_qss_quadratic_coupling(reactions_variables, mech_variables, qss_variables, exit_error=exit_error)

    return coupled_qss_matrix


def needs_dictionaries(qss_variables, coupling_variables):
    """Computes the need dependencies between qss species

    Parameters
    ----------
    qss_variables : dict
        dictionary of QSS related parameters needed for computation
    coupling_variables : dict
        dictionary of variables related to the coupling of QSS species

    Returns
    -------
    needs : dict
        which other species are needed by one
    needs_full : dict
        the uncompacted version of the prior, meaning there is no group present
    is_needed : dict
        by which a species is needed
    needs_count : dict
        how many species does one need

    """
    species_qss_names = qss_variables['species_qss_names']
    qss_index_local = qss_variables['qss_index_local']

    coupled_qss_matrix = qss_variables['coupled_qss_matrix']
    coupled_qss_matrix_modified = coupled_qss_matrix.copy()

    group = coupling_variables['group']

    # Setting the lines of species belonging to the same group to be equal

    if group:
        for key in group:
            sum_line = 0
            sum_column = 0
            for value in group[key]:
                sum_line += coupled_qss_matrix_modified[qss_index_local[value], :]
                sum_column += coupled_qss_matrix_modified[:, qss_index_local[value]]
            for value in group[key]:
                coupled_qss_matrix_modified[qss_index_local[value], :] = sum_line
                coupled_qss_matrix_modified[:, qss_index_local[value]] = sum_column
                for value2 in group[key]:
                    coupled_qss_matrix_modified[qss_index_local[value], qss_index_local[value2]] = 0
                    coupled_qss_matrix_modified[qss_index_local[value2], qss_index_local[value]] = 0
    coupled_qss_matrix_modified[coupled_qss_matrix_modified >= 1] = 1

    index_spec_line, index_spec_column = np.nonzero(coupled_qss_matrix)

    needs = {}
    needs_full = {}
    needs_count = {}
    is_needed = {}

    # needs dictionary
    for index_spec_needs, spec_needs in enumerate(species_qss_names):
        token_needs_full = []
        for index_spec_is_needed in index_spec_column[index_spec_line == index_spec_needs]:
            token_needs_full.append(species_qss_names[index_spec_is_needed])
        needs_full[spec_needs] = token_needs_full

    index_spec_line, index_spec_column = np.nonzero(coupled_qss_matrix_modified)

    for index_spec_needs, spec_needs in enumerate(species_qss_names):
        token_needs = []
        for index_spec_is_needed in index_spec_column[index_spec_line == index_spec_needs]:
            not_in_a_group = True
            if group:
                for key in group:
                    if species_qss_names[index_spec_is_needed] in group[key]:
                        not_in_a_group = False
                        if key not in token_needs:
                            token_needs.append(key)

            if not_in_a_group and species_qss_names[index_spec_is_needed] not in token_needs:
                token_needs.append(species_qss_names[index_spec_is_needed])

        not_in_a_group = True
        if group:
            for key, values in list(group.items()):
                if spec_needs in values:
                    needs[key] = token_needs
                    needs_count[key] = len(token_needs)
                    not_in_a_group = False

        if not_in_a_group:
            needs[spec_needs] = token_needs
            needs_count[spec_needs] = len(token_needs)

    # is_needed dictionary
    for index_spec_is_needed, spec_is_needed in enumerate(species_qss_names):
        token_is_needed = []
        for index_spec_needs in index_spec_line[index_spec_column == index_spec_is_needed]:
            not_in_a_group = True
            if group:
                for key in group:
                    if species_qss_names[index_spec_needs] in group[key]:
                        not_in_a_group = False
                        if key not in token_is_needed:
                            token_is_needed.append(key)

            if not_in_a_group and species_qss_names[index_spec_needs] not in token_is_needed:
                token_is_needed.append(species_qss_names[index_spec_needs])

        not_in_a_group = True
        if group:
            for key, values in list(group.items()):
                if spec_is_needed in values:
                    is_needed[key] = token_is_needed
                    not_in_a_group = False

        if not_in_a_group:
            is_needed[spec_is_needed] = token_is_needed

    # Checking groups
    species_in_groups = []
    which_group = {}
    for key in group:
        species_in_groups += group[key]
        for spec in group[key]:
            which_group[spec] = key

    species_not_in_groups = [spec for spec in species_qss_names if spec not in species_in_groups]

    # If a species needs a group and is needed by it, it is added to the group
    for ite in range(3):
        for spec in species_qss_names:
            if spec not in species_in_groups:
                for key in group:
                    if key in needs[spec] and spec in needs[key]:

                        # Updating group species
                        group[key].append(spec)

                        # Updating needs
                        needs[key].remove(spec)
                        needs.pop(spec)

                        # Updating is_needed
                        is_needed[key].remove(spec)
                        is_needed.pop(spec)

                        # Updating needs_count
                        needs_count.pop(spec)
                        needs_count[key] -= 1

                        # Updating species sorting
                        species_in_groups.append(spec)
                        species_not_in_groups.remove(spec)

                        # Changing it in other species data by the group
                        for second_spec in species_not_in_groups:
                            if spec in needs[second_spec]:
                                needs[second_spec].remove(spec)
                                needs[second_spec].append(key)
                                if second_spec not in is_needed[key]:
                                    is_needed[key].append(second_spec)
                            if spec in is_needed[second_spec]:
                                is_needed[second_spec].remove(spec)
                                is_needed[second_spec].append(key)
                                if second_spec not in needs[key]:
                                    needs[key].append(second_spec)
                                    needs_count[key] += 1
                        break

        # Checking for groups to merge
        changing_size = True
        while changing_size:
            copy_group = group.copy()
            changing_size = False
            for key in copy_group:
                if key in group:
                    group_needs = needs[key]
                    for item in group_needs:
                        if item.startswith('group'):
                            if key in needs[item]:

                                # Merging the groups
                                group[key] += group[item]
                                group.pop(item)

                                # Updating needs
                                needs[key].remove(item)
                                needs.pop(item)

                                # Updating is_needed
                                is_needed[key].remove(item)
                                is_needed.pop(item)

                                # Updating needs_count
                                needs_count.pop(item)
                                needs_count[key] -= 1

                                changing_size = True

    return needs, needs_full, is_needed, needs_count


def coupled_species_group(qss_variables, maintain_given_order=True):
    """Groups species that are coupled together

    Parameters
    ----------
    qss_variables : dict
        dictionary of QSS related parameters needed for computation
    maintain_given_order : bool
        if True, the species inside the groups will be sorted according to the global list (Default value = True)

    Returns
    -------
    group : dict
        dictionary of grouped species

    """

    species_qss_names = qss_variables['species_qss_names']
    coupled_qss_matrix = qss_variables['coupled_qss_matrix']

    index_spec_line, index_spec_column = np.nonzero(coupled_qss_matrix)

    group = {}

    count = 0

    for index_spec_needs in range(len(coupled_qss_matrix[0])):
        buffer_list = []
        for index_spec_is_needed in index_spec_column[index_spec_line == index_spec_needs]:
            # if i needs j and j needs i --> group
            if coupled_qss_matrix[index_spec_needs, index_spec_is_needed] == 1 \
                    and coupled_qss_matrix[index_spec_is_needed, index_spec_needs] == 1:
                buffer_list.append(species_qss_names[index_spec_needs])
                buffer_list.append(species_qss_names[index_spec_is_needed])

                # for index_indirect_spec_needs in index_spec_column[index_spec_line == index_spec_is_needed]:
                #     # if another species k is needed by i and needs j or needs i and is needed by j --> group
                #     if (coupled_qss_matrix[index_indirect_spec_needs, index_spec_is_needed] == 1
                #         and coupled_qss_matrix[index_spec_is_needed, index_spec_needs] == 1) \
                #             or (coupled_qss_matrix[index_spec_needs, index_spec_is_needed] == 1
                #                 and coupled_qss_matrix[index_spec_is_needed, index_indirect_spec_needs] == 1):
                #         buffer_list.append(species_qss_names[index_indirect_spec_needs])

                group['group_' + str(count)] = buffer_list
                count += 1

    # Cleaning groups
    for key in group:
        group[key] = list(set(group[key]))

    # Merging groups with common species (2 loops to ensure that there is no duplicate)
    for loop in range(2):
        group_temp = {}
        group_left = list(group.keys())
        for key in group:
            new_group_values = []
            group_taken = []
            for key2 in group_left:
                intersection = np.intersect1d(group[key], group[key2])
                if len(intersection) > 0:
                    new_group_values += group[key2]
                    group_taken.append(key2)
            if new_group_values:
                group_taken.append(key)
                for key in group_taken:
                    if key in group_left:
                        group_left.remove(key)

                group_temp[key] = list(set(new_group_values))

        group = group_temp

    # Cleaning empty groups
    for key in list(group.keys()):
        for key2 in list(group.keys()):
            if key != key2 and key in list(group.keys()) and set(np.intersect1d(group[key], group[key2])) == set(
                    group[key]):
                del group[key]

    # Renumbering
    buffer_group = {}
    for index, key in enumerate(group):
        buffer_group['group_' + str(index)] = group[key]
    group = buffer_group

    # Reordering species inside group
    if maintain_given_order:
        for key in group:
            buffer_list = []
            for spec in species_qss_names:
                if spec in group[key]:
                    buffer_list.append(spec)
            group[key] = buffer_list

    return group


def reorder(coupling_variables):
    """Function reordering the QSS dependencies matrix into a sub-triangular matrix

    Parameters
    ----------
    coupling_variables : dict
        dictionary of variables related to the coupling of QSS species

    Returns
    -------
    species_index : int
        index of the species in the reordered order

    """
    needs_count = coupling_variables['needs_count']
    is_needed = coupling_variables['is_needed']

    decouple_index = {}

    k = 0
    needs_count_regress = needs_count.copy()

    counter = 0
    while len(decouple_index) < len(needs_count):
        counter += 1
        if counter > len(needs_count) + 1:
            break
        needs_count_base = needs_count_regress.copy()
        for i in needs_count_base:
            if needs_count_regress[i] == 0:
                decouple_index[i] = k  # index for the qss reordering
                del needs_count_regress[i]  # delete the species from the needs dictionary
                for j in is_needed[i]:
                    if j in list(needs_count_regress.keys()):
                        needs_count_regress[j] -= 1  # lower the number of needs for the species who needed species i
                k += 1

    decouple_index_keys = list(decouple_index.keys())
    decouple_index_values = list(decouple_index.values())
    species_index = dict(list(zip(decouple_index_values, decouple_index_keys)))

    if len(species_index) != len(needs_count):
        logger.warning('Huge cluster of coupled QSS species')
        species_index = {}

    return species_index


def qss_coeffs(print_variables, constants, mech_variables, reactions_variables, qss_variables, coupling_variables,
               spec_qss):
    """Prints the expressions of the coefficients needed for the computation of QSS species concentrations

    Parameters
    ----------
    print_variables : dict
        dictionary containing the f90 file to write in and code specific infos
    constants : dict
        dictionary of constants needed for computation
    mech_variables : dict
        dictionary of mechanism related parameters needed for computation
    reactions_variables : dict
        dictionary of reactions related parameters needed for computation
    qss_variables : dict
        dictionary of qss related parameters needed for computation
    coupling_variables : dict
        dictionary of variables related to the coupling of qss species
    spec_qss : str
        QSS species for which the coefficients must be written


    """
    f = print_variables['f']
    precision = print_variables['precision']

    ns = constants['ns']

    mech = mech_variables['mech']
    nur = mech_variables['nur']
    nup = mech_variables['nup']

    reac_direction = reactions_variables['reac_direction']
    reac_type = reactions_variables['reac_type']
    is_reversible = reactions_variables['is_reversible']
    reverse_index = reactions_variables['reverse_index']

    species_qss_names = qss_variables['species_qss_names']
    qss_index_global = qss_variables['qss_index_global']
    qss_index_local = qss_variables['qss_index_local']
    index_coupled_qss_i = qss_variables['index_coupled_qss_i']
    index_coupled_qss_j = qss_variables['index_coupled_qss_j']
    coupled_qss_matrix = qss_variables['coupled_qss_matrix']

    group = coupling_variables['group']
    species_in_group = coupling_variables['species_in_group']
    which_group = coupling_variables['which_group']

    deltaReactantsR = abs(nur)
    deltaProductsR = abs(nup)

    j = 0
    for i in range(ns):
        if i not in qss_index_global:
            deltaReactantsR = np.delete(deltaReactantsR, j, 0)
            deltaProductsR = np.delete(deltaProductsR, j, 0)
        else:
            j += 1

    for i in range(deltaReactantsR.shape[0]):
        for j in range(deltaReactantsR.shape[1]):
            if deltaReactantsR[i, j] == - deltaProductsR[i, j]:
                deltaReactantsR[i, j] = 0
                deltaProductsR[i, j] = 0

    deltaReactantsR_qss = deltaReactantsR
    deltaProductsR_qss = deltaProductsR

    indrri, indrrj = np.nonzero(deltaReactantsR_qss)
    indrpi, indrpj = np.nonzero(deltaProductsR_qss)

    text = ''

    i = qss_index_local[spec_qss]

    # Denominator
    if indrrj[indrri == i].any() or indrpj[indrpi == i].any():
        i = int(i)
        text += '    ' + species_qss_names[i] + '_denom = tiny(1.0_' + precision + ') + ( 0.0_' + precision + ' '
        # Forward
        for j in indrrj[indrri == i]:
            reactants = mech.reaction(j).reactants
            reactants_keys = list(reactants.keys())
            reactants_keys = tools.convert_to_valid_fortran_var(reactants_keys)
            reactants_values = list(reactants.values())
            text += '&' + '\n' + '             + k(r' + str(j + 1) + reac_direction[j] + ')'
            for k in range(len(reactants_keys)):
                if reactants_keys[k] != spec_qss:
                    c_name = 'c'
                    s_name = 's'
                    if reactants_keys[k] in species_qss_names:
                        c_name = 'cqss'
                        s_name = 'sqss'
                        if reactants_keys[k] in species_in_group:
                            if spec_qss in species_in_group \
                                    and which_group[reactants_keys[k]] != which_group[spec_qss]:
                                if reactants_values[k] != 1:
                                    text += ' *' + c_name + '(' + s_name + \
                                            reactants_keys[k] + ')** ' + str(reactants_values[k]) + '_' + precision
                                else:
                                    text += '* ' + c_name + '(' + s_name + reactants_keys[k] + ') '
                        else:
                            if reactants_values[k] != 1:
                                text += ' *' + c_name + '(' + s_name + \
                                        reactants_keys[k] + ')** ' + str(reactants_values[k]) + '_' + precision
                            else:
                                text += '* ' + c_name + '(' + s_name + reactants_keys[k] + ') '
                    else:
                        if reactants_values[k] != 1:
                            text += ' *' + c_name + '(' + s_name + \
                                    reactants_keys[k] + ')** ' + str(reactants_values[k]) + '_' + precision
                        else:
                            text += '* ' + c_name + '(' + s_name + reactants_keys[k] + ') '
            if reac_type[j] == 2:
                text += '* m(mM' + str(j + 1) + ') '
        # Backward
        for j in indrpj[indrpi == i]:
            if is_reversible[j] == 1:
                reactants = mech.reaction(j).products
                reactants_keys = list(reactants.keys())
                reactants_keys = tools.convert_to_valid_fortran_var(reactants_keys)
                reactants_values = list(reactants.values())
                text += '&' + '\n' + '             + k(r' + str(j + 1) + reac_direction[reverse_index[j]] + ')'
                for k in range(len(reactants_keys)):
                    if reactants_keys[k] != spec_qss:
                        c_name = 'c'
                        s_name = 's'
                        if reactants_keys[k] in species_qss_names:
                            c_name = 'cqss'
                            s_name = 'sqss'
                            if reactants_keys[k] in species_in_group:
                                if spec_qss in species_in_group \
                                        and which_group[reactants_keys[k]] != which_group[spec_qss]:
                                    if reactants_values[k] != 1:
                                        text += ' *' + c_name + '(' + s_name + \
                                                reactants_keys[k] + ')** ' + str(reactants_values[k]) + '_' + precision
                                    else:
                                        text += '* ' + c_name + '(' + s_name + reactants_keys[k] + ') '
                            else:
                                if reactants_values[k] != 1:
                                    text += ' *' + c_name + '(' + s_name + \
                                            reactants_keys[k] + ')** ' + str(reactants_values[k]) + '_' + precision
                                else:
                                    text += '* ' + c_name + '(' + s_name + reactants_keys[k] + ') '
                        else:
                            if reactants_values[k] != 1:
                                text += ' *' + c_name + '(' + s_name + \
                                        reactants_keys[k] + ')** ' + str(reactants_values[k]) + '_' + precision
                            else:
                                text += '* ' + c_name + '(' + s_name + reactants_keys[k] + ') '
                if reac_type[j] == 2:
                    text += '* m(mM' + str(j + 1) + ') '

        text += '  )\n\n'

        # Numerator
        text += '    ' + species_qss_names[i] + '_num = ( 0.0_' + precision + ' '
        # Forward
        for j in indrpj[indrpi == i]:
            reactants = mech.reaction(j).reactants
            reactants_keys = list(reactants.keys())
            reactants_keys = tools.convert_to_valid_fortran_var(reactants_keys)
            reactants_values = list(reactants.values())

            products = mech.reaction(j).products
            products_keys = list(products.keys())
            products_keys = tools.convert_to_valid_fortran_var(products_keys)
            products_values = list(products.values())
            products = dict(list(zip(products_keys, products_values)))

            spec_of_group = []
            if spec_qss in species_in_group:
                spec_of_group = group[which_group[spec_qss]]

            if not set(reactants_keys).intersection(spec_of_group):
                if products[species_qss_names[i]] != 1:
                    text += '&' + '\n' + '     +' + str(
                            products[species_qss_names[i]]) + '_' + precision + ' * k(r' + str(
                            j + 1) + \
                            reac_direction[j] + ')'
                else:
                    text += '&' + '\n' + '             + k(r' + str(j + 1) + reac_direction[j] + ')'
                for k in range(len(reactants_keys)):
                    if reactants_keys[k] != spec_qss:
                        c_name = 'c'
                        s_name = 's'
                        if reactants_keys[k] in species_qss_names:
                            c_name = 'cqss'
                            s_name = 'sqss'
                        if reactants_values[k] != 1:
                            text += ' *' + c_name + '(' + s_name + \
                                    reactants_keys[k] + ')** ' + str(reactants_values[k]) + '_' + precision
                        else:
                            text += '* ' + c_name + '(' + s_name + reactants_keys[k] + ') '
                if reac_type[j] == 2:
                    text += '* m(mM' + str(j + 1) + ') '
        # Backward
        for j in indrrj[indrri == i]:
            if is_reversible[j] == 1:
                reactants = mech.reaction(j).products
                reactants_keys = list(reactants.keys())
                reactants_keys = tools.convert_to_valid_fortran_var(reactants_keys)
                reactants_values = list(reactants.values())

                products = mech.reaction(j).reactants
                products_keys = list(products.keys())
                products_keys = tools.convert_to_valid_fortran_var(products_keys)
                products_values = list(products.values())
                products = dict(list(zip(products_keys, products_values)))
                spec_of_group = []
                if spec_qss in species_in_group:
                    spec_of_group = group[which_group[spec_qss]]

                if not set(reactants_keys).intersection(spec_of_group):
                    if products[species_qss_names[i]] != 1:
                        text += '&' + '\n' + '     +' + str(
                                products[species_qss_names[i]]) + '_' + precision + ' * k(r' + str(j + 1) + \
                                reac_direction[reverse_index[j]] + ')'
                    else:
                        text += '&' + '\n' + '             + k(r' + str(j + 1) + reac_direction[reverse_index[j]] + ')'
                    for k in range(len(reactants_keys)):
                        if reactants_keys[k] != spec_qss:
                            c_name = 'c'
                            s_name = 's'
                            if reactants_keys[k] in species_qss_names:
                                c_name = 'cqss'
                                s_name = 'sqss'
                            if reactants_values[k] != 1:
                                text += ' *' + c_name + '(' + s_name + \
                                        reactants_keys[k] + ')** ' + str(reactants_values[k]) + '_' + precision
                            else:
                                text += '* ' + c_name + '(' + s_name + reactants_keys[k] + ') '
                    if reac_type[j] == 2:
                        text += '* m(mM' + str(j + 1) + ') '

        text += '  )\n\n'

        text += '    ' + species_qss_names[i] + '_ct = ' + species_qss_names[i] + '_num / ' + species_qss_names[
            i] + '_denom\n\n'

        # Coupled QSS species coefficient

        for j in index_coupled_qss_j[index_coupled_qss_i == i]:
            if (coupled_qss_matrix[i, j] == 1 and coupled_qss_matrix[j, i] == 1) or \
                    (coupled_qss_matrix[i, j] == 1 and spec_qss in species_in_group and species_qss_names[
                        j] in species_in_group):
                text += '    ' + species_qss_names[i] + '_' + species_qss_names[j] + ' = - ( 0.0_' + precision + ' '
                # Forward
                for k in indrpj[indrpi == i]:
                    reactants = mech.reaction(k).reactants
                    reactants_keys = list(reactants.keys())
                    reactants_keys = tools.convert_to_valid_fortran_var(reactants_keys)
                    reactants_values = list(reactants.values())

                    products = mech.reaction(k).products
                    products_keys = list(products.keys())
                    products_keys = tools.convert_to_valid_fortran_var(products_keys)
                    products_values = list(products.values())
                    products = dict(zip(products_keys, products_values))

                    if species_qss_names[j] in reactants_keys:
                        text += '&' + '\n' + '          + k(r' + str(k + 1) + reac_direction[k] + ') '
                        if products[species_qss_names[i]] != 1:
                            text += ' * ' + str(products[species_qss_names[i]]) + '_' + precision + ' '

                        for s in range(len(reactants_keys)):
                            if reactants_keys[s] != spec_qss:
                                c_name = 'c'
                                s_name = 's'
                                if reactants_keys[s] in species_qss_names:
                                    c_name = 'cqss'
                                    s_name = 'sqss'
                                    if reactants_keys[s] in species_in_group:
                                        if spec_qss in species_in_group and which_group[reactants_keys[s]] != \
                                                which_group[spec_qss]:
                                            if reactants_values[s] != 1:
                                                text += '*' + c_name + '(' + s_name + \
                                                        reactants_keys[s] + ')' + '** ' + str(reactants_values[
                                                                                                  s]) + '_' + precision
                                            else:
                                                text += '* ' + c_name + '(' + s_name + reactants_keys[s] + ') '
                                    else:
                                        if reactants_values[s] != 1:
                                            text += '*' + c_name + '(' + s_name + \
                                                    reactants_keys[s] + ')** ' + str(
                                                    reactants_values[s]) + '_' + precision
                                        else:
                                            text += '* ' + c_name + '(' + s_name + reactants_keys[s] + ') '
                                else:
                                    if reactants_values[s] != 1:
                                        text += '*' + c_name + '(' + s_name + \
                                                reactants_keys[s] + ')** ' + str(reactants_values[s]) + '_' + precision
                                    else:
                                        text += '* ' + c_name + '(' + s_name + reactants_keys[s] + ') '
                        if reac_type[k] == 2:
                            text += '* m(mM' + str(k + 1) + ') '

                # Backward
                for k in indrrj[indrri == i]:
                    if is_reversible[k] == 1:
                        reactants = mech.reaction(k).products
                        reactants_keys = list(reactants.keys())
                        reactants_keys = tools.convert_to_valid_fortran_var(reactants_keys)
                        reactants_values = list(reactants.values())

                        products = mech.reaction(k).reactants
                        products_keys = list(products.keys())
                        products_keys = tools.convert_to_valid_fortran_var(products_keys)
                        products_values = list(products.values())
                        products = dict(zip(products_keys, products_values))

                        if species_qss_names[j] in reactants_keys:
                            text += '&' + '\n' + '          + k(r' + str(k + 1) + reac_direction[
                                reverse_index[k]] + ') '
                            if products[species_qss_names[i]] != 1:
                                text += ' * ' + str(products[species_qss_names[i]]) + '_' + precision + ' '

                            for s in range(len(reactants_keys)):
                                if reactants_keys[s] != spec_qss:
                                    c_name = 'c'
                                    s_name = 's'
                                    if reactants_keys[s] in species_qss_names:
                                        c_name = 'cqss'
                                        s_name = 'sqss'
                                        if reactants_keys[s] in species_in_group:
                                            if spec_qss in species_in_group and which_group[reactants_keys[s]] != \
                                                    which_group[spec_qss]:
                                                if reactants_values[s] != 1:
                                                    text += '*' + c_name + '(' + s_name + \
                                                            reactants_keys[s] + ')** ' + str(
                                                            reactants_values[s]) + '_' + precision
                                                else:
                                                    text += '* ' + c_name + '(' + s_name + reactants_keys[s] + ') '
                                        else:
                                            if reactants_values[s] != 1:
                                                text += '*' + c_name + '(' + s_name + \
                                                        reactants_keys[s] + ')** ' + str(
                                                        reactants_values[s]) + '_' + precision
                                            else:
                                                text += '* ' + c_name + '(' + s_name + reactants_keys[s] + ') '
                                    else:
                                        if reactants_values[s] != 1:
                                            text += '*' + c_name + '(' + s_name + \
                                                    reactants_keys[s] + ')** ' + str(
                                                    reactants_values[s]) + '_' + precision
                                        else:
                                            text += '* ' + c_name + '(' + s_name + reactants_keys[s] + ') '
                            if reac_type[k] == 2:
                                text += '* m(mM' + str(k + 1) + ') '

                text += ' ) / ' + species_qss_names[i] + '_denom\n\n'

    else:
        text += '    ' + species_qss_names[i] + '_ct = 0.0_' + precision + '\n\n'

    f.write(text)


def solve_qss_concentration(print_variables, constants, mech_variables, reactions_variables, qss_variables,
                            coupling_variables):
    """Prints the expressions of QSS species concentrations computation solutions

    Parameters
    ----------
    print_variables : dict
        dictionary containing the f90 file to write in and code specific infos
    constants : dict
        dictionary of constants needed for computation
    mech_variables : dict
        dictionary of mechanism related parameters needed for computation
    reactions_variables : dict
        dictionary of reactions related parameters needed for computation
    qss_variables : dict
        dictionary of QSS related parameters needed for computation
    coupling_variables : dict
        dictionary of variables related to the coupling of QSS species

    """

    f = print_variables['f']
    precision = print_variables['precision']

    species_qss_names = qss_variables['species_qss_names']
    qss_index_local = qss_variables['qss_index_local']
    coupled_qss_matrix = qss_variables['coupled_qss_matrix']

    group = coupling_variables['group']
    species_index = coupling_variables['species_index']
    concentration = {}
    species_qss_ordered = []

    # Error case
    if not species_index:
        group = {'group_0': species_qss_names}
        species_index = {0: 'group_0'}

    V_name = []
    C_name = []
    V_text = {}
    C_text = {}

    count = 0
    for n in range(len(species_index)):
        i = species_index[n]
        if i in list(group.keys()):

            s = group[i]

            Variable = [['0.0_' + precision] * len(s) for i in range(len(s))]
            Concentration = ['0.0_' + precision] * len(s)
            Constant = ['0.0_' + precision] * len(s)

            for l in range(len(s)):
                l = int(l)
                Constant[l] = s[l] + '_ct'
                Concentration[l] = 'cqss(sqss' + s[l] + ')'
                for m in range(len(s)):
                    m = int(m)
                    if l == m:
                        Variable[l][m] = '1.0_' + precision
                    else:
                        if coupled_qss_matrix[qss_index_local[s[l]], qss_index_local[s[m]]] == 1:
                            Variable[l][m] = s[l] + '_' + s[m]

            V, X, C, names = qss_expressions_gauss_pivoting(print_variables, Variable, Concentration, Constant,
                                                            qss_group_number=count)
            count += 1
            V_name += names['A_name']
            C_name += names['B_name']

            V_text[i] = names['A_val']
            C_text[i] = names['B_val']

            for k in range(len(s)):
                concentration[s[len(s) - 1 - k]] = X[len(s) - 1 - k]
                species_qss_ordered.append(s[len(s) - 1 - k])
        else:
            concentration[i] = i + '_ct'
            species_qss_ordered.append(i)

    species_in_group = []  # list of species belonging to groups
    which_group = {}  # dictionary of which group correspond to a species
    for g in group:
        for spec in group[g]:
            which_group[spec] = g
            if spec not in species_in_group:
                species_in_group.append(spec)

    coupling_variables['species_in_group'] = species_in_group
    coupling_variables['which_group'] = which_group

    already_written = []
    for variable in V_name + C_name:
        f.write('    real(pr) :: ' + variable + '\n')

    f.write('\n')

    used_groups = []

    for index_spec, speqss in enumerate(species_qss_ordered):
        count = 0
        if speqss not in already_written:
            if speqss in species_in_group:
                for spec in group[which_group[speqss]]:
                    qss_coeffs(print_variables, constants, mech_variables, reactions_variables, qss_variables,
                               coupling_variables, species_qss_ordered[index_spec + count])
                    already_written.append(species_qss_ordered[index_spec + count])
                    count += 1
                count = 0

                if which_group[speqss] not in used_groups:
                    for values in V_text[which_group[speqss]]:
                        f.write('    ' + values + '\n')

                    for values in C_text[which_group[speqss]]:
                        f.write('    ' + values + '\n')

                    used_groups.append(which_group[speqss])

                f.write('\n')

                for spec in group[which_group[speqss]]:
                    text = '    cqss(sqss' + species_qss_ordered[index_spec + count] + ') = ' + concentration[
                        species_qss_ordered[index_spec + count]]
                    text += '  \n\n'
                    count += 1
                    f.write(text)
            else:
                qss_coeffs(print_variables, constants, mech_variables, reactions_variables, qss_variables,
                           coupling_variables, speqss)
                text = '    cqss(sqss' + speqss + ') = ' + concentration[speqss]
                text += '  \n\n'
                f.write(text)

    text = """\
    cqss = max(cqss, tiny(1.0_pr))
    cqss = min(cqss, 1e03_pr)

    ! Necessary for very low concentration variations
    ! Needs to be integrated only for Cantera
    do index=1,nqss
      if (cqss(index) /= cqss(index)) then
        cqss(index) = tiny(1.0_pr)
      end if
    end do

    return
  end subroutine get_qss


"""

    f.write(text)


def print_qss_concentrations(print_variables, constants, mech_variables, reactions_variables, qss_variables,
                             force=True, maintain_given_order=True):
    """Function calling all other functions needed for the QSS handling

    Parameters
    ----------
    print_variables : dict
        dictionary containing the f90 file to write in and code specific infos
    constants : dict
        dictionary of constants needed for computation
    mech_variables : dict
        dictionary of mechanism related parameters needed for computation
    reactions_variables : dict
        dictionary of reactions related parameters needed for computation
    qss_variables : dict
        dictionary of reactions related parameters needed for computation
    force : bool
        if True, species and reactions will be removed from the mechanism if they induce QSS error (Default value = True)
    maintain_given_order : bool
        if True, the species inside the groups will be sorted according to the global list (Default value = True)

    """
    coupled_qss_matrix = qss_coupling(constants, mech_variables, reactions_variables, qss_variables, force=force)

    qss_variables['coupled_qss_matrix'] = coupled_qss_matrix

    group = coupled_species_group(qss_variables, maintain_given_order=maintain_given_order)

    coupling_variables = {}
    coupling_variables['group'] = group

    needs, needs_full, is_needed, needs_count = needs_dictionaries(qss_variables, coupling_variables)

    coupling_variables['needs'] = needs
    coupling_variables['needs_full'] = needs_full
    coupling_variables['needs_count'] = needs_count
    coupling_variables['is_needed'] = is_needed

    print_constants_qss(print_variables, constants, qss_variables)

    species_index = reorder(coupling_variables)

    coupling_variables['species_index'] = species_index

    solve_qss_concentration(print_variables, constants, mech_variables, reactions_variables, qss_variables,
                            coupling_variables)


def new_variable_name_generator(rank=0, position=0, iteration=0, upper_case=False):
    """Creates dynamically a variable name

    Parameters
    ----------
    rank : int
        number of the group to be processed (Default value = 0)
    position : int
        position in the linearized matrix (Default value = 0)
    iteration : int
        iteration of the Gauss elimination (Default value = 0)
    upper_case : bool
        if True, chooses the upper cases (Default value = False)

    Returns
    -------
    name : str
        generated variable name

    """

    name_parts = []
    indices = [rank + 1, position, iteration + 1]

    for index in indices:

        if upper_case:
            letters = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
        else:
            letters = list('abcdefghijklmnopqrstuvwxyz')

        if index > 26:
            div = index // 26
            rem = index % 26
            if div > 26:
                div_div = div // 26
                rem_div = div % 26
                new_letter = letters[div_div - 1] + letters[rem_div - 1] + letters[rem - 1]
            else:
                new_letter = letters[div - 1] + letters[rem - 1]
            letter = new_letter
        else:
            letter = letters[index - 1]

        name_parts.append(letter)

    name = name_parts[0] + '_' + name_parts[1] + '_' + name_parts[2]

    return name


def qss_expressions_gauss_pivoting(print_variables, A, X0, B, qss_group_number=0):
    r"""Solves the system :math:`\mathbf{A}\mathbf{x_0} = \mathbf{b}` by a Gauss elimination and returns
    :math:`\mathbf{A}`, :math:`\mathbf{x}` and :math:`\mathbf{b}`.

    Note
    ----

        Notations:

        - :math:`QSS` is the QSS species
        - :math:`P_i` and :math:`R_i` are the non-QSS species associated to the :math:`i^\text{th}`
          reaction :math:`(R_i)`. The former are associated to the **reactants** while the latter are the **products**
        - :math:`\mathbf{A}` is the matrix of the **coupling coefficients between QSS species**
        - :math:`\mathbf{x}` is the array of the **species concentrations**
        - :math:`\mathbf{b}` is the array of the **constant** parts


    Example
    -------
    The following example deals with 5 reaction with 3 QSS species. Let :math:`QSS_i` be the three QSS species,
    for :math:`i \in \{1,2,3\}`.

    Given by these set of reactions:

    .. math::
        (R_1): QSS_3 + R_1 \quad \stackrel{k_1^{f,b}}{\Leftrightarrow} &\quad  QSS_1 + P_1 \\\\
        (R_2): QSS_2 + R_2  \quad \stackrel{k_2^f}{\Rightarrow} & \quad QSS_1 + P_2 \\\\
        (R_3): QSS_3 + R_3  \quad \stackrel{k_3^f}{\Rightarrow} & \quad QSS_2 + P_3 \\\\
        (R_4): R_4 \quad \stackrel{k_4^f}{\Rightarrow} & \quad QSS_1 + P_4 \\\\
        (R_5): R_5 \quad \stackrel{k_5^f}{\Rightarrow} & \quad QSS_3 + P_5

    The system is expressed to identify as :math:`\mathbf{A}\mathbf{x} = \mathbf{b} \ (^*)`:

    .. math::
        \begin{equation*}
            (^*) =
            \begin{pmatrix}
                \frac{\mathrm{d}[QSS_1]}{\mathrm{d}t} = k_1^f[QSS_3][*] - k_1^b[QSS_1][*] + k_2^f[QSS_2][*] + k_4^f[*] \\\\
                \frac{\mathrm{d}[QSS_2]}{\mathrm{d}t} = - k_2^f[QSS_2][*] + k_3^f[QSS_3][*] \\\\
                \frac{\mathrm{d}[QSS_3]}{\mathrm{d}t} = - k_1^f[QSS_3][*] + k_1^b[QSS_1][*] - k_3^f[QSS_3][*] + k_5^f[*]
            \end{pmatrix}
        \end{equation*}

    According to QSS hypothesis, :math:`\frac{\mathrm{d}[QSS_i]}{\mathrm{d}t} = 0`, for :math:`\in \{1,2,3\}`,
    therefore :math:`(^*)` becomes:

    .. math::
        \begin{equation*}
            \mathbf{A}\mathbf{x} = \mathbf{b} \Leftrightarrow
            \begin{pmatrix}
                [QSS_1] - \frac{k_2^f[*]}{k_1^b[*]} [QSS_2] - \frac{k_1^f[*]}{k_1^b[*]} [QSS_3]
                    =   \frac{k_4^f[*]}{k_1^b[*]} \\\\
                [QSS_2] - \frac{k_3^f[*]}{k_2^f[*]} [QSS_3]
                    =   0 \\\\
                - \frac{k_1^b[*]}{k_1^f[*] + k_3^f[*]} [QSS_1] + [QSS_3]
                    =   \frac{k_5^f[*]}{k_1^f[*] + k_3^f[*]}
            \end{pmatrix}
        \end{equation*}

    with,

    .. math::
        \begin{equation*}
            \mathbf{A} =
            \begin{pmatrix}
                1               & f(QSS_1,QSS_2)    & f(QSS_1,QSS_3)  \\\\
                0               & 1                 & f(QSS_2,QSS_3)  \\\\
                f(QSS_1,QSS_3)  & 0                 & 1
            \end{pmatrix}
        \end{equation*}

    .. math::
        \begin{equation*}
            \mathbf{x} =
            \begin{pmatrix}
                [QSS_1]\\\\
                [QSS_2]\\\\
                [QSS_3]
            \end{pmatrix}
        \end{equation*}

    .. math::
        \begin{equation*}
            \mathbf{b} =
            \begin{pmatrix}
                QSS_1^{ct}\\\\
                QSS_2^{ct}\\\\
                QSS_3^{ct}
            \end{pmatrix}
        \end{equation*}


    Parameters
    ----------
    print_variables :
        dictionary containing the f90 file to write in and code specific infos
    A : list
        :math:`n\times n` list of strings
    X0 : list
        :math:`n` list with the variables names
    B : list
        :math:`n` list of strings
    qss_group_number : int
        QSS group number (Default value = 0)

    Returns
    -------
    A : list
        modified expression of matrix :math:`\mathbf{A}` after applying Gauss pivot
    X : list
        :math:`\mathbf{x}`
    B : list
        modified expression of matrix :math:`\mathbf{B}` after applying Gauss pivot
    name : str

    """

    precision = print_variables['precision']

    names = {}
    A_name = []
    B_name = []
    A_val = []
    B_val = []

    # System dimension
    n = len(A[0])

    X = [''] * n
    for i in range(n):
        X[i] = 'X' + str(i)

        delta_A = np.zeros([n, n])
    for i in range(n):
        for j in range(n):
            if A[i][j] != '0.0_' + precision:
                delta_A[i, j] = 1

    for k in range(n - 1):  # loop through the diagonal except last value
        pivot = A[k][k]
        if pivot == '0.0_' + precision:  # if the pivot is nul, switches current line with the following
            temp = np.array(A[k + 1][:])
            A[k + 1][:] = A[k][:]
            A[k][:] = temp

            temp = str(B[k + 1])
            B[k + 1] = B[k]
            B[k] = temp

            pivot = A[k][k]

        for i in range(k, n - 1):  # loop through the lines following the pivot line
            num = A[i + 1][k]  # numerator to have the element under the pivot equal to 0
            B = list(B)
            if num != '0.0_' + precision:
                if pivot != '1.0_' + precision:
                    if num != '1.0_' + precision:
                        B[i + 1] = ' (' + B[i + 1] + ') - (' + B[
                            int(k)] + ') * (' + num + ') / (' + pivot + ')'
                else:
                    if num != '1.0_' + precision:
                        B[i + 1] = ' (' + B[i + 1] + ') - (' + B[int(k)] + ') * (' + num + ')'
                    else:
                        B[i + 1] = ' (' + B[i + 1] + ') - (' + B[int(k)] + ')'

            var_name = new_variable_name_generator(qss_group_number, i + 1, k, upper_case=True)
            B_val.append(var_name + ' = ' + B[i + 1])
            B_name.append(var_name)
            B[i + 1] = var_name

            # Set matrix column of the previous pivot to zero
            if k > 0:
                delta_A[:, k - 1] = 0
            indi, indj = np.nonzero(delta_A)

            for j in indj[indi == k]:
                if A[i + 1][j] != '0.0_' + precision:
                    if num != '0.0_' + precision:
                        if pivot != '1.0_' + precision:
                            if num != '1.0_' + precision:
                                A[i + 1][j] = ' (' + A[i + 1][j] + ') - (' \
                                              + A[k][j] + ') * (' + num + ') / (' + pivot + ')'
                        else:
                            if num != '1.0_' + precision:
                                A[i + 1][j] = ' (' + A[i + 1][j] + ') - (' + A[k][
                                    j] + ') * (' + num + ')'
                            else:
                                A[i + 1][j] = ' (' + A[i + 1][j] + ') - (' + A[k][j] + ')'
                else:
                    if num != '0.0_' + precision:
                        if pivot != '1.0_' + precision:
                            if num != '1.0_' + precision:
                                A[i + 1][j] = ' - (' \
                                              + A[k][j] + ') * (' + num + ') / (' + pivot + ')'
                        else:
                            if num != '1.0_' + precision:
                                A[i + 1][j] = ' - (' + A[k][j] + ') * (' + num + ')'
                            else:
                                A[i + 1][j] = ' - (' + A[k][j] + ')'

                if j != k and A[i + 1][j] != '0.0_' + precision:
                    var_name = new_variable_name_generator(qss_group_number, ((i + 1) * n) + j, k)
                    A_name.append(var_name)
                    A_val.append(var_name + ' = ' + A[i + 1][j])
                    A[i + 1][j] = var_name

            A[i + 1][k] = '0.0_' + precision

    # System resolution giving the values of the X array

    for i in range(n):
        X = list(X)
        B[i] = str(B[i])

    n = n - 1
    if A[n][n] != '1.0_' + precision:
        X[n] = '( ' + B[n] + ' ) / ( ' + A[n][n] + ' )'
    else:
        X[n] = B[n]
    for i in range(1, n + 1):
        sumprod = ''
        for j in range(i):
            if A[n - i][n - j] != '0.0_' + precision:
                if A[n - i][n - j] == '1.0_' + precision:
                    sumprod += ' - ' + str(X0[n - j])
                else:
                    sumprod += ' - (' + A[n - i][n - j] + ') * ' + X0[n - j]

        if sumprod == '':
            if A[n - i][n - i] != '1.0_' + precision:
                X[n - i] = B[n - i] + ' / ' + A[n - i][n - i]
            else:
                X[n - i] = B[n - i]
        else:
            if A[n - i][n - i] == '1.0_' + precision:
                X[n - i] = B[n - i] + sumprod
            else:
                X[n - i] = '(' + B[n - i] + sumprod + ') / (' + A[n - i][n - i] + ')'

    # Lines length correction
    for index_x, x in enumerate(X):
        size_x = len(x)
        index = 0
        threshold = 40
        while index + threshold < size_x:
            if x[index + threshold:].find(') ') > -1:
                index_space = x[index + threshold:].find(') ') + index + threshold + 1
                x = x[:index_space] + ' &\n                ' + x[index_space:]
                index = index_space + 21
                size_x = len(x)
            else:
                break
        X[index_x] = x

    names['A_name'] = A_name
    names['B_name'] = B_name
    names['A_val'] = A_val
    names['B_val'] = B_val

    return A, X, B, names


def check_qss_quadratic_coupling(reactions_variables, mech_variables, qss_variables, exit_error=False):
    """Function to find quadratic coupling between QSS species and propose a better subset of QSS species

    Parameters
    ----------
    reactions_variables : dict
        dictionary of reactions related parameters needed for computation
    mech_variables : dict
        dictionary of mechanism related parameters needed for computation
    qss_variables : dict
        dictionary of qss related parameters needed for computation
    exit_error : bool
        if True exits the python in case of error (Default value = False)

    Returns
    -------
    problematic_reactions : list
        list of problematic reactions
    species_qss_names_proposition :
        list of quadratic coupling-free QSS species

    """

    nur = mech_variables['nur']
    nup = mech_variables['nup']

    is_reversible = reactions_variables['is_reversible']

    species_qss_names = qss_variables['species_qss_names']
    qss_index_global = qss_variables['qss_index_global']
    index_coupled_qss_i = qss_variables['index_coupled_qss_i']
    index_coupled_qss_j = qss_variables['index_coupled_qss_j']
    index_qss_reactions_i = qss_variables['index_qss_reactions_i']
    index_qss_reactions_j = qss_variables['index_qss_reactions_j']

    species_qss_names_proposition = species_qss_names.copy()
    problematic_species = {}
    problematic_reactions = []
    fail = 0

    for i in range(len(species_qss_names)):
        for j in index_coupled_qss_j[index_coupled_qss_i == i]:
            if j != i:
                for k in index_qss_reactions_j[index_qss_reactions_i == i]:
                    if (nur[qss_index_global[i], k] >= 1 and nur[qss_index_global[j], k] >= 1) or \
                            ((nup[qss_index_global[i], k] >= 1 and nup[qss_index_global[j], k] >= 1) and
                             is_reversible[k] == 1):
                        if k not in problematic_reactions:
                            problematic_reactions.append(k)
                            logger.debug(
                                '!!! WARNING !!! No quadratic coupling allowed (species ' + species_qss_names[i]
                                + ' and ' + 'species ' + species_qss_names[j] + ' in reaction ' + str(k + 1) + ')')
                        # Dictionnary on how many times a species is involved in a conflict
                        if not species_qss_names_proposition[i] in problematic_species:
                            problematic_species[species_qss_names[i]] = 0
                        else:
                            problematic_species[species_qss_names[i]] += 1

                        if not species_qss_names_proposition[j] in problematic_species:
                            problematic_species[species_qss_names[j]] = 0
                        else:
                            problematic_species[species_qss_names[j]] += 1
                        fail = 1

        # Species coupled with themselves
        for k in index_qss_reactions_j[index_qss_reactions_i == i]:
            if nur[qss_index_global[i], k] > 1 or (nup[qss_index_global[i], k] > 1 and is_reversible[k] == 1):
                if k not in problematic_reactions:
                    problematic_reactions.append(k)
                    logger.debug('!!! WARNING !!! No quadratic coupling allowed (species ' + species_qss_names[i]
                                 + ' with coefficient ' + str(int(nur[qss_index_global[i], k])) + ' in reaction ' + str(
                            k + 1) + ')')
                # Dictionnary on how many times a species is involved in a conflict
                if not species_qss_names_proposition[i] in problematic_species:
                    problematic_species[species_qss_names[i]] = 0
                else:
                    problematic_species[species_qss_names[i]] += 1

                fail = 1

    # Proposition of a better QSS set
    while problematic_species:
        species_qss_names_proposition.remove(
                max(problematic_species.keys(), key=(lambda key: problematic_species[key])))
        problematic_species = {}

        for i in range(len(species_qss_names)):
            for j in index_coupled_qss_j[index_coupled_qss_i == i]:
                if j != i \
                        and species_qss_names[i] in species_qss_names_proposition \
                        and species_qss_names[j] in species_qss_names_proposition:
                    for k in index_qss_reactions_j[index_qss_reactions_i == i]:
                        if (nur[qss_index_global[i], k] >= 1 and nur[qss_index_global[j], k] >= 1) or \
                                ((nup[qss_index_global[i], k] >= 1 and nup[qss_index_global[j], k] >= 1) and
                                 is_reversible[k] == 1):
                            # Dictionary on how many times a species is involved in a conflict
                            if not species_qss_names[i] in problematic_species:
                                problematic_species[species_qss_names[i]] = 0
                            else:
                                problematic_species[species_qss_names[i]] += 1

                            if not species_qss_names[j] in problematic_species:
                                problematic_species[species_qss_names[i]] = 0
                            else:
                                problematic_species[species_qss_names[j]] += 1

    if exit_error and fail == 1:
        logger.error('Maybe you should try this quadratic coupling-free subset:' + str(species_qss_names_proposition))
        exit('!!! ERROR !!! The set of QSS species provided is not valid ! Check previous messages for intel')

    return problematic_reactions, species_qss_names_proposition
