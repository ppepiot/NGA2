"""Module implementing DRGEP for species and reactions reduction"""

import ARCANE.sampling as sampling
import ARCANE.display as display
import ARCANE.kwdict as kwdict
import ARCANE.analysis as analysis
import ARCANE.postproc as postproc

import numpy as np

import sys


logger = display.Logger()
logger.set_log('logReduction')
kwdict = kwdict.Kwdict()


def drgep_ranked(cases_list, mechanism, reduction_type, integrity=True, sensitivity=0):
    """Sort species using DRGEP method.DRGEP is applied to all samples in sdb, uses mechanism in mechanism.


    Parameters
    ----------
    cases_list :
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    reduction_type :
        type of reduction (species/S or reactions/R)
    integrity :
        check if there is no mass sink or source (Default value = True)
    sensitivity :
        Activate aiding DRGEP with SA (Sensitivity analysis)
        (0 to be not activated, a value to be the sensor on the DRGEP coefficients. See sensitivity_analysis in analysis.) (Default value = 0)

    Returns
    -------
    array containing ranked species
        
        
    Created: 17/11/14 [PP]

    Last modified: 19/08/30 [JW]

    """

    if reduction_type not in ['species', 'S', 'reactions', 'R']:
        logger.error(" The type of reduction must be either 'species' (or 'S') or  'reactions' (or 'R')")
        sys.exit()

    ctmech = mechanism.ctmech
    species_names = ctmech.species_names
    reactions = ctmech.reactions()

    myns = mechanism.ns
    mynr = mechanism.nr

    if reduction_type in ['species', 'S']:
        EP = np.zeros(myns, 'd')
    else:
        EP = np.zeros(mynr, 'd')

    # Putting together cases with the same targets
    cases_dict = {}
    for case in cases_list:
        targets = repr(case.targets)
        if targets not in list(cases_dict.keys()):
            cases_dict[targets] = [case]
        else:
            cases_dict[targets].append(case)

    cases_tuples = list(zip(cases_dict.keys(), cases_dict.values()))

    for targets, cases in cases_tuples:

        targets = eval(targets)

        # Create sample database on cases
        sampled_dict = sampling.create_samples_database(cases, mechanism, ['T', 'P', 'Y', 'HR'])

        EPtmp = error_propagation(mechanism, sampled_dict, targets, reduction_type)
        EP = np.fmax(EP, EPtmp)

    if integrity and reduction_type in ['species', 'S']:
        no_more_consumed, no_more_produced, remove_if_removed = integrity_check(mechanism)
        # Integrity check
        if integrity:
            for spec_removed in no_more_produced:
                for spec_impacted in no_more_produced[spec_removed]:
                    EP[ctmech.species_index(spec_impacted)] = EP[ctmech.species_index(spec_removed)]

            for spec_removed in no_more_consumed:
                for spec_impacted in no_more_consumed[spec_removed]:
                    max_EP_spec_impacted = [EP[ctmech.species_index(spec_impacted)] \
                                            for spec_impacted in no_more_consumed[spec_removed]]
                    max_value = np.max([EP[ctmech.species_index(spec_removed)],
                                        np.max(max_EP_spec_impacted)])
                    EP[ctmech.species_index(spec_removed)] = max_value
                    EP[ctmech.species_index(spec_impacted)] = max_value

    if reduction_type in ['species', 'S']:
        # Keep important species
        important_species = important_species_identification(cases_list, mechanism)

        # Add sensitivity to the important species
        # Species are switched according to an arbitrary coefficient defined in sensitivity
        # => If the species is a target, the DRGEP coefficient will be 1.
        # => If the species is under the sensitivity coefficient, treatment is done to keep the DRGEP in this order.
        # => If the species is above the sensitivity coefficient, treatment is done to put the species coefficient
        # bigger than the DRGEP.
        if sensitivity:
            logger.info('Sensitivity applied for species coefficient above ' + str(sensitivity))
            for id_case, case in enumerate(cases_list):
                result_analysis = analysis.sensitivity_analysis(case, reduction_type)[case.mechanism.name][case.myid]
                for id_err, err in enumerate(list(case.error_dict.keys())):
                    if id_case + id_err == 0:
                        sa_species = result_analysis[err]
                    else:
                        sa_species = {key_dic: max(value_dic, result_analysis[err][key_dic])\
                                      for key_dic, value_dic in sa_species.items()}
            for id_spec, spec in enumerate(species_names):
                if spec in important_species:
                    EP[ctmech.species_index(spec)] = 1
                elif sa_species[id_spec] < sensitivity:
                    EP[ctmech.species_index(spec)] = EP[ctmech.species_index(spec)]*max(min(sa_species.values()), 1e-30)\
                                                     /max(sa_species.values())
                else:
                    EP[ctmech.species_index(spec)] = sa_species[id_spec]/max(sa_species.values())

                #logger.info('Sensible species :' + ', '.join([sa_species[spec] for spec in species_names if sa_species[spec] < sensitivity]))
        # Classical method
        else:
            EP = [1.0 if spec in important_species else EP[ctmech.species_index(spec)] for spec in species_names]

        drgep_dict = dict(zip(species_names, EP))

    elif reduction_type in ['reactions', 'R']:

        # Add sensitivity to the important reactions
        # Species are switched according to an arbitrary coefficient defined in sensitivity
        # => If the reaction is a target, the DRGEP coefficient will be 1.
        # => If the reaction is under the sensitivity coefficient, treatment is done to keep the DRGEP in this order.
        # => If the reaction is above the sensitivity coefficient, treatment is done to put the species coefficient
        # bigger than the DRGEP.

        if sensitivity:
            logger.info('Sensitivity applied for reactions coefficient above ' + str(sensitivity))
            for id_case, case in enumerate(cases_list):
                result_analysis = analysis.sensitivity_analysis(case, reduction_type)[case.mechanism.name][case.myid]
                for id_err, err in enumerate(list(case.error_dict.keys())):
                    if id_case + id_err == 0:
                        sa_reactions = result_analysis[err]
                    else:
                        sa_reactions = {key_dic: max(value_dic, result_analysis[err][key_dic])\
                                      for key_dic, value_dic in sa_reactions.items()}
            for id_reac, reac in enumerate(reactions):
                if sa_reactions[reac.equation] < sensitivity:
                    EP[id_reac] = EP[id_reac]*max(min(sa_reactions.values()), 1e-30)\
                                                     /max(sa_reactions.values())
                else:
                    EP[id_reac] = sa_reactions[reac.equation]/max(sa_reactions.values())

        drgep_dict = dict(zip(reactions, EP))

    else:
        drgep_dict = {}
        logger.error('DRGEP is not applicable on this case !')

    return drgep_dict


def error_propagation(mechanism, sampled_dict, targets, reduction_type):
    """Computes error propagation coefficients.

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    sdb :
        sample database on which DRGEP is applied
    targets :
        list of reduction targets
    reduction_type :
        type of reduction performed (species or reactions)

    Returns
    -------
    dictionary of species or reactions with their corresponding coefficient


    Created: 17/11/14 [PP]

    Last modified: 18/04/03 [QC]

    """

    myns = mechanism.ns
    mynr = mechanism.nr

    if reduction_type in ['species', 'S']:
        EP = np.zeros(myns, 'd')
    elif reduction_type in ['reactions', 'R']:
        EP = np.zeros(mynr, 'd')
    else:
        EP = []
        logger.error('Wrong reduction type !')

    data_dict = postproc.extend_data_dict(sampled_dict)
    DIC_spec, DIC_reac = compute_DIC(mechanism, sampled_dict, reduction_type)

    alpha_norm = scaling_coefficient(mechanism, data_dict, targets)

    logger.terminator('\r')
    n_data = len(sampled_dict['grid'])

    for i_sample in range(n_data):
        logger.info('####### Samples treated : ' + str(i_sample + 1) + ' / ' + str(n_data) + ' #######')

        local_DIC_spec = DIC_spec[:, :, i_sample]
        if reduction_type in ['reactions', 'R']:
            local_DIC_reac = DIC_reac[:, :, i_sample]
        else:
            local_DIC_reac = None

        EPtmp = local_error_propagation(mechanism, local_DIC_spec, local_DIC_reac, targets,
                                        reduction_type, alpha_norm[:, i_sample])

        EP = np.fmax(EP, EPtmp)
    logger.info('\n')
    logger.terminator('\n')

    return EP


def local_error_propagation(mechanism, local_DIC_spec, local_DIC_reac, targets, reduction_type, alpha_norm_loc):
    """Computes error propagation coefficients on one sample.

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    mysample :
        sample on which error propagation is applied
    targets :
        list of reduction targets
    reduction_type :
        type of reduction performed (species or reactions)
    alpha_norm_loc :
        scaling coefficient for the given sample
   
    Returns
    -------
    list of error coefficients

    Created: 17/11/14 [PP]

    Last modified: 19/01/28 [QC]

    """

    myns = mechanism.ns
    mynr = mechanism.nr
    net = mechanism.network
    ctmech = mechanism.ctmech
    species_names = ctmech.species_names

    # Parameters
    EPmin = 1e-7
    EPcomp = np.zeros(myns, 'd')
    EP_spec = np.zeros(myns, 'd')
    EP_reac = np.zeros(mynr, 'd')
    EP_spec_dict = {}

    # Indices of targets
    target_indices = [species_names.index(target) for target in targets 
                            if target not in ['HeatRelease', 'HeatRelease']]
    if 'HeatRelease' in targets or 'HeatRelease' in targets:
        target_indices.append(myns)

    # ------------
    # Error Propagation (EP)
    # ------------

    # Go through each target
    for index_target_local, index_target_global in enumerate(target_indices):

        # Initialize working array - Works for HR as well
        EPtmp = local_DIC_spec[index_target_global, :] * alpha_norm_loc[index_target_local]

        # Initialize EP arrays
        array_up = np.zeros(myns, 'd')
        array_down = np.zeros(myns, 'd')

        # Initialize array_up
        array_up[:] = -1.0
        for i in range(myns):
            if i == index_target_global or EPtmp[i] < EPmin:
                continue

            array_up[i] = EPtmp[i]

        # Iterate until all relevant coefficients have been included

        flag = True
        while flag:

            # Init inner loop
            flag = False
            array_down[:] = -1.0

            # Loop over array_up
            for i in range(myns):
                indj = net.indsj[net.indsi == i]
                # If coeff is positive, store coeff and continue

                if array_up[i] > 0.0:
                    coeff_up = array_up[i]

                    # Loop over all species
                    for j in indj:
                        coeff_down = local_DIC_spec[i, j] * coeff_up

                        # New coeff counts if i!=j and > EPmin
                        if i != j and coeff_down > EPmin:
                            flag = True
                            # Update EPtmp and array_down for next iteration
                            if coeff_down > EPtmp[j]:
                                EPtmp[j] = coeff_down
                                array_down[j] = coeff_down

            if list(EPcomp) == list(EPtmp):
                flag = False
            else:
                EPcomp = EPtmp

            if targets[index_target_local] != 'HeatRelease':
                EP_spec_dict[species_names[index_target_global]] = EPtmp
            else:
                EP_spec_dict['HeatRelease'] = EPtmp

            array_up[:] = array_down[:]
        # Adjust maximum coefficient for target/species pair
        EP_spec = np.fmax(EP_spec, EPtmp)

    if reduction_type in ['reactions', 'R']:
        EPtmp = np.zeros(mynr, 'd')

        for index_target_local, index_target_global in enumerate(target_indices):
            # Initialize working array - Works for HR as well
            R_target = np.zeros(mynr, 'd')
            for index_spec in range(myns):
                if targets[index_target_local] == 'HeatRelease':
                    Rtmp = EP_spec_dict['HeatRelease'][index_spec] * local_DIC_reac[index_spec, :]
                else:
                    Rtmp = EP_spec_dict[species_names[index_target_global]][index_spec] \
                           * local_DIC_reac[index_spec, :]

                R_target = np.fmax(R_target, Rtmp)

            EPtmp[:] = R_target[:]

            # Adjust maximum coefficient for target/species pair
            EP_reac = np.fmax(EP_reac, EPtmp)

    if reduction_type in ['species', 'S']:
        EP = EP_spec
    else:
        EP = EP_reac

    return EP

def compute_DIC(mechanism, data_dict, reduction_type):
    """Computes dictionary.

    Parameters
    ----------
    mechanism :
        param data_dict
    data_dict :
        dictionary of data to be converted
    reduction_type :
        type of reduction performed (species or reactions)


    """

    # Parameters
    n_data = len(data_dict['grid'])
    n_spec = mechanism.ns
    n_reac = mechanism.nr
    net = mechanism.network
    nup = net.nup
    nur = net.nur
    nu = nup - nur

    # Summed values for normalisation
    HR_prod_reac = np.array([sum(data_dict['hr reactions'][:, index_data][data_dict['hr reactions'][:, index_data] > 0])
                    for index_data in range(n_data)])
    HR_cons_reac = np.array([sum(data_dict['hr reactions'][:, index_data][data_dict['hr reactions'][:, index_data] < 0])
                    for index_data in range(n_data)])

    index_HR = n_spec

    # Workspace
    DIC_spec = np.zeros([n_spec + 1, n_spec, n_data], 'd')

    # Evaluate DIC(i,j)
    for i in range(n_spec):
        # reactions containing species i
        booli = net.indr[net.inds == i]
        indj = net.indsj[net.indsi == i]

        for j in indj:
            boolj = net.indr[net.inds == j]  # reactions containing species j
            indk = np.intersect1d(booli, boolj)  # reactions containing species i and j

            # Compute the DIC
            for k in indk:
                DIC_spec[i, j, :] += (nup[i, k] - nur[i, k]) * data_dict['net rates of progress'][k, :]

            # Normalize
            DIC_spec[i, j, :] = abs(DIC_spec[i, j, :]) / np.maximum(np.maximum(data_dict['creation rates'][i, :],
                                                              data_dict['destruction rates'][i, :]), 1e-60)

        # Heat release term
        DIC_spec[index_HR, i, :] += abs(data_dict['hr species'][i, :])

        # Normalize
        DIC_spec[index_HR, i, :] = abs(DIC_spec[index_HR, i, :]) \
                                   / np.maximum(np.maximum(abs(HR_prod_reac), abs(HR_cons_reac)), 1e-60)

    if reduction_type in ['reactions', 'R']:
        # Workspace
        DIC_reac = np.zeros([n_spec + 1, n_reac, n_data], 'd')

        # Evaluate DIC(i,j)
        for i in range(n_spec):
            # reactions containing species i
            indk = net.indr[net.inds == i]  # reactions containing species i

            # Compute the DIC
            for k in indk:
                DIC_reac[i, k, :] = (nup[i, k] - nur[i, k]) * data_dict['net rates of progress'][k, :]

                # Normalize
                DIC_reac[i, k, :] = abs(DIC_reac[i, k, :]) / np.maximum(np.maximum(data_dict['creation rates'][i, :],
                                                              data_dict['destruction rates'][i, :]), 1e-60)

        for i in range(n_reac):

            # Heat release term
            DIC_reac[index_HR, i] += abs(data_dict['hr reactions'][i, :])

            # Normalize
            DIC_reac[index_HR, i] = abs(DIC_reac[index_HR, i]) \
                                   / np.maximum(np.maximum(abs(HR_prod_reac), abs(HR_cons_reac)), 1e-60)

    else:
        DIC_reac = None

    return DIC_spec, DIC_reac


def integrity_check(mechanism):
    """Dictionaries with intel about integrity.
    
    Gives dictionaries telling:
    - which species are no more consumed if a given one (the key of the dict) is removed
    - which species are no more produced if a given one (the key of the dict) is removed
    - which species it would be wise to remove if a given one (the key of the dict) is removed

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
        
    Returns
    -------
    return the 3 dictionaries previously described in the same order

    Created: 18/04/03 [QC]
    """
    myns = mechanism.ns
    net = mechanism.network
    ctmech = mechanism.ctmech
    species_names = ctmech.species_names
    nup = net.nup
    nur = net.nur
    nu = nup - nur

    # Integrity dict: grouping species that induce truncated path
    no_more_produced = {}
    no_more_consumed = {}
    # Global dictionary giving the species that needs to be removed if key is removed
    remove_if_removed = {}

    for i in range(myns):
        nu_integrity = nu.copy()
        # reactions containing species i
        indj = net.indsj[net.indsi == i]

        for j in indj:

            # Integrity check
            if [x for x in nu_integrity[j, :] if x != 0]:
                if np.max([x for x in nu_integrity[j, :] if
                           x != 0]) < 0:  # if species i is removed species j is no more produced
                    if species_names[i] not in no_more_produced:
                        no_more_produced[species_names[i]] = []
                    else:
                        if species_names[j] not in no_more_produced[species_names[i]]:
                            no_more_produced[species_names[i]].append(species_names[j])
                    # Global dict
                    if species_names[i] not in remove_if_removed:
                        remove_if_removed[species_names[i]] = []
                    else:
                        if species_names[j] not in remove_if_removed[species_names[i]]:
                            remove_if_removed[species_names[i]].append(species_names[j])

                elif np.min([x for x in nu_integrity[j, :] if
                             x != 0]) > 0:  # if species i is removed species j is no more consumed
                    if species_names[i] not in no_more_consumed:
                        no_more_consumed[species_names[i]] = []
                    else:
                        if species_names[j] not in no_more_consumed[species_names[i]]:
                            no_more_consumed[species_names[i]].append(species_names[j])
                    # Global dict
                    if species_names[i] not in remove_if_removed:
                        remove_if_removed[species_names[i]] = []
                    else:
                        if species_names[j] not in remove_if_removed[species_names[i]]:
                            remove_if_removed[species_names[i]].append(species_names[j])

    return no_more_consumed, no_more_produced, remove_if_removed


# Needs to be made more efficient
def scaling_coefficient(mechanism, samples_dict, targets):
    """
    Computes scaling coefficient for DRGEP.

    Parameters
    ----------
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object
    samples_dict :
        sample database on which DRGEP is applied
    targets :
        list of reduction targets

    Returns
    -------
    number_of_targets times number_of_samples matrix of coefficients

    
    Created: 18/04/03 [QC]

    """

    ctmech = mechanism.ctmech
    species_names = ctmech.species_names

    elements = ctmech.element_names

    # Indices of targets
    target_indices = [species_names.index(target) for target in targets if target not in ['HeatRelease', 'HeatRelease']]
    target_indices.append(len(species_names))

    # Heat release
    try:
        i_HR = targets.index('HeatRelease')
    except ValueError:
        i_HR = None

    # Store composition of species in a n_elements*n_species matrix
    element_in_species = np.zeros([len(elements), len(species_names)])
    for index_spec, spec in enumerate(species_names):
        comp = ctmech.species(index_spec).composition
        for index_element, element in enumerate(elements):
            if element in comp:
                element_in_species[index_element, index_spec] = comp[element]

    # Initialize workspace
    n_samples = len(samples_dict['grid'])
    n_elements = len(elements)
    n_targets = len(targets)

    # Identifying the different cases in the samples
    cases = list(set([samples_dict['case id'][index] for index in range(n_samples)]))
    n_cases = len(cases)

    ind_e, ind_s = np.nonzero(element_in_species)
    alpha = np.zeros([n_elements, n_targets, n_samples])
    alpha_temp = np.zeros([n_elements, n_targets, n_samples])
    max_alpha = np.ones([n_elements, n_targets, n_cases])  # Ones instead of zeros to avoid zero division
    alpha_norm = np.zeros([n_targets, n_samples])

    for i_sample in range(n_samples):

        # Sample case index
        case_index = cases.index(samples_dict['case id'][i_sample])

        # ------------
        # Scaling coefficient
        # ------------

        # Element production
        for index_element, element in enumerate(elements):
            P_element = 0.0
            # Species containing the element
            spec_with_element = ind_s[ind_e == index_element]

            for index_spec in spec_with_element:
                # P_e = Sum(N_e,s * max(P_s - C_s, 0))
                P_element += max(element_in_species[index_element, index_spec] * (
                        samples_dict['destruction rates'][index_spec, i_sample]
                        - samples_dict['creation rates'][index_spec, i_sample]), 1e-60)

            for index_target, target in enumerate(targets):
                if target != 'HeatRelease':
                    # alpha_e,t = (N_e,s * |P_t - C_t|) / P_e
                    alpha[index_element, index_target, i_sample] = element_in_species[index_element, int(
                        target_indices[index_target])] * abs(samples_dict['destruction rates'][target_indices[index_target], i_sample]
                            - samples_dict['creation rates'][target_indices[index_target], i_sample]) / P_element

                    if alpha[index_element, index_target, i_sample] > max_alpha[index_element, index_target, case_index]:

                        max_alpha[index_element, index_target, case_index] = alpha[index_element, index_target, i_sample]

        if i_HR is not None:
            alpha[0, i_HR, i_sample] = abs(samples_dict['HR'][i_sample])
            if alpha[0, i_HR, i_sample] > max_alpha[0, i_HR, case_index]:
                max_alpha[0, i_HR, case_index] = alpha[0, i_HR, i_sample]

    # Reconstructing
    max_alpha_temp = np.zeros([n_elements, n_targets, n_samples])
    for i_sample in range(n_samples):
        # Sample case index
        case_index = cases.index(samples_dict['case id'][i_sample])
        max_alpha_temp[:, :, i_sample] = max_alpha[:, :, case_index]

    max_alpha = max_alpha_temp

    for index_element, element in enumerate(elements):
        for index_target, target in enumerate(targets):
            alpha_temp[index_element, index_target, :] = alpha[index_element, index_target, :] \
                                                             / max_alpha[index_element, index_target, :]

    # Normalization: alpha_t = max_e(alpha_e,t/max_sample(alpha_e,t))
    for index_target, target in enumerate(targets):
        alpha_norm[index_target, :] = np.max(alpha_temp[:, index_target, :], 0)

    return alpha_norm


def important_species_identification(cases_list, mechanism):
    """
    Identifies which the important species (fuel, oxidizer, products and diluents).
    Those species are the one with a mass faction representing more than 1% of the fresh and burnt gases.

    Parameters
    ----------
    cases_list :
        list of class :func:`~ARCANE.cases.Case` objects
    mechanism :
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    list of important species

    """

    species_names = mechanism.ctmech.species_names

    important_species = []

    for case in cases_list:
        data_dict = case.data_dict(mechanism)
        for spec in species_names:
            if spec not in important_species:
                y_start = data_dict[spec][0]
                y_end = data_dict[spec][-1]

                if y_start > 1e-3 or y_end > 2e-2:
                    important_species.append(spec)

                if spec in case.targets:
                    important_species.append(spec)

    return important_species

