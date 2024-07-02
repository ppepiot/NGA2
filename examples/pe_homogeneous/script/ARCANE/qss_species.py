"""Module implementing LOI for species reduction"""

import ARCANE.sampling as sampling
import ARCANE.drgep as drgep
import ARCANE.postproc as postproc
import ARCANE.tools as tools
import ARCANE.display as display

import numpy as np

logger = display.Logger()
logger.set_log('logReduction')


def qss_ranked_species_LOI(cases_list, mechanism):
    """Sorts species using LOI method.LOI is applied to all samples in sdb, uses mechanism in mechanism

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects on which LOI is applied
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    array containing ranked species


    Created: 18/01/22 [QC]
        
    Last modified: 18/03/13 [QC]

    """

    myns = mechanism.ns
    ctmech = mechanism.ctmech
    species_names = ctmech.species_names

    LOI = np.zeros(myns, 'd')
    LOI_temp = np.zeros(myns, 'd')

    # Putting together cases with the same targets
    cases_dict = {}
    for case in cases_list:
        targets = repr(case.targets)
        if targets not in cases_dict:
            cases_dict[targets] = [case]
        else:
            cases_dict[targets].append(case)

    cases_tuples = list(zip(cases_dict.keys(), cases_dict.values()))

    for targets, cases in cases_tuples:

        targets = eval(targets)

        # Create sample database on cases
        samples_dict = sampling.create_samples_database(cases, mechanism)
        n_data = len(samples_dict['grid'])

        samples_dict = postproc.extend_data_dict(samples_dict)
        DIC_spec, DIC_reac = drgep.compute_DIC(mechanism, samples_dict, 'S')
        alpha_norm = drgep.scaling_coefficient(mechanism, samples_dict, targets)

        logger.terminator('\r')
        for i_sample in range(n_data):
            logger.info('####### Samples treated : ' + str(i_sample + 1) + ' / ' + str(n_data) + ' #######')

            local_DIC_spec = DIC_spec[:, :, i_sample]
            local_DIC_reac = None

            EP = drgep.local_error_propagation(mechanism, local_DIC_spec, local_DIC_reac, targets,
                                        'S', alpha_norm[:, i_sample])

            local_dict = {}
            local_dict['T'] = samples_dict['T'][i_sample]
            local_dict['P'] = samples_dict['P'][i_sample]
            local_dict['Y'] = samples_dict['Y'][i_sample]

            timescale = tools.local_timescale(local_dict, mechanism, concentration_weighting=True)

            LOItmp = np.multiply(np.sqrt(EP), timescale)
            LOI_temp = np.fmax(LOI_temp, LOItmp)

        LOI = np.fmax(LOI, LOI_temp)

    # Keep important species
    important_species = drgep.important_species_identification(cases_list, mechanism)
    LOI = [1.0 if spec in important_species else LOI[ctmech.species_index(spec)] for spec in species_names]

    LOI_dict = dict(zip(ctmech.species_names, LOI))

    logger.terminator('\n')

    return LOI_dict


def qss_ranked_species_AOI(cases_list, mechanism):
    """Sorts species using max of mass fractions integral "Area Of Importance (AOI)".
    Gets the corresponding dumped scalars and sort it according to their maximum.

    Parameters
    ----------
    cases_list : list
        list of class :func:`~ARCANE.cases.Case` objects on which LOI is applied
    mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
        class :func:`~ARCANE.mechanisms.Mechanism` object

    Returns
    -------
    array containing ranked species
        

    Created: 18/01/22 [QC]
        
        
    Last modified: 18/03/13 [QC]

    """

    ctmech = mechanism.ctmech
    species = ctmech.species_names

    # Create sample database on cases
    scalars_db = sampling.scalars_database(cases_list, mechanism)

    max_integral_dict = {}
    for i, spec in enumerate(species):
        max_integral_dict[spec] = np.max(scalars_db[:, i])

    AOI = {}
    for i, spec in enumerate(species):
        AOI[spec] = max_integral_dict[spec]

    important_species = drgep.important_species_identification(cases_list, mechanism)

    for i, spec in enumerate(species):
        if AOI[spec] == 0 or spec in important_species:
            AOI[spec] = 1.0

    return AOI
