"""Carrying the full reduction of a mechanism

Author: QC (2018/06/18)
"""

import os
import time

import numpy as np

import ARCANE.cases as cases
import ARCANE.custom.custom_kinetics as custom
import ARCANE.display as display
import ARCANE.drgep as drgep
import ARCANE.error as error
import ARCANE.graphs as graphs
import ARCANE.kwdict as kwdict
import ARCANE.lumping as lumping
import ARCANE.mechanisms as mechanisms
import ARCANE.qss_species as qss
import ARCANE.tools as tools
import ARCANE.database as database

from matplotlib.pyplot import cm
import pandas as pd
from graphviz import Source
import matplotlib

logger = display.Logger()
logger.set_log('logReduction')

kwdict = kwdict.Kwdict()
graph = graphs.Graph()

class Automatic:

    def __init__(self, cases_list, mechanism):
        """Class for automatic reduction of a mechanism on a set of cases

        Parameters
        ----------
        cases_list : list
            list of class :func:`~ARCANE.cases.Case` objects
        mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            class :func:`~ARCANE.mechanisms.Mechanism` object

        """

        # List of cases for the reduction
        self.cases_list = cases_list
        # Mechanism object to be reduced
        self.mechanism = mechanism
        # Mechanism object to be compared to
        self.super_root_mechanism = mechanism
        # Flag for remote installation
        self.remote_install = True
        # Flag for overwriting cases profiles
        self.overwrite = False
        # 'list' or 'mechanism' depending on where the non valid QSS should be removed
        self.remove_from = 'mechanism'
        # Diluent that won't be discarded even if not used
        self.diluent = ['N2']
        # Threshold value for stopping the reduction
        self.error_threshold = 1
        # Threshold relative to the error or absolute
        self.relative_threshold = True
        # Flag for using AOI identification of QSS
        self.fast_qss_identification = False
        # Method for case restoration
        self.restore = 'auto'
        # Plotting lists
        self.xdata = []
        self.ydata = []
        # Reduction directory
        self.reddir = 'reduced'
        # Takes 'first' or 'last' valid mechanism
        self.which_reduced = 'last'
        # Error factors for each reduction type
        self.error_factors = [0, 0, 0, 0]
        # Step for discarding species or reactions
        self.step = "auto"
        # Flag for skipping non valid reduction steps
        # Can be a list or a boolean
        self.non_stop = [True, False, True]
        # Flag for automatically removing all null values in DRGEP
        self.skip_zeros = False
        # Flag for going through ranked items with unitary value
        self.bypass_ones = False
        # Flag for plotting the reduction graphs
        self.reduction_graph = 0
        # Flag for checking if compensation
        self.compensation = 0
        # Flag for dichotomy applied
        self.dichotom = False
        # Criterion on the DRGEP or LOI/AOI values that are considered too high
        # Allows to save some computation time (only activated if non_stop is too)
        self.stop_criterion = True
        self.criterion_factors = [1, 0.1, 10]
        # Initializes the counter for the species and the reactions
        self.counter_spec = 0
        self.counter_reac = 0
        # Default plotting values
        self.error_plot = True
        self.max_error_plot = True
        # Pressure taken for pdep reactions simplification
        self.pdep_value = 1e5
        # List of tuples for the automatic stepping  of species
        self.ns_auto_stepping = [(300, 30), (200, 20), (100, 10), (50, 1)]
        # List of tuples for the automatic stepping  of reactions
        self.nr_auto_stepping = [(10000, 500), (5000, 200), (1000, 50), (500, 10), (200, 1)]
        # Order of the reverse Arrhenius fit and it parameters
        self.fit_order = 0
        self.start_temperature = 300
        self.end_temperature = 2500
        self.number_of_temperatures = 6000
        self.fit_pressure = 101325
        # Testing the fortran conversion of the skeletal
        self.testing_reverse = True

        # List to store errors
        self.max_errors = []
        self.max_error_spec = []
        self.max_error_reac = []
        self.max_error_qss = []

        # Variables used for reduction
        self.relevant_var = ''
        self.name_of_keys = ''
        self.name_of_method = ''
        self.type_index = -1
        self.ranked_dict = {}

        # reference for timing the reduction
        self.start_time = time.time()
        self.cantera_run_time = 0

    def set_mechanism(self, mechanism):
        """Sets the ARCANE mechanism that will be used for reduction

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        """

        self.mechanism = mechanism

    def set_which_reduced(self, position):
        """Sets if the mechanism picked as reduced is the one before the error goes above tolerance
        or the last one satisfying this tolerance

        Parameters
        ----------
        position : str
            Either 'first' or 'last'. The former will returns the first mechanism going above tolerance
            while the latter returns the last one satisfying this tolerance

        """

        self.which_reduced = position

    def set_super_root_mechanism(self, mechanism):
        """Sets the super root ARCANE mechanism that will be used for reduction. This mechanism is independent
        of the reduction step. This will avoid error addition when performing serial reductions

        Parameters
        ----------
        mechanism :
            class :func:`~ARCANE.mechanisms.Mechanism` object

        """

        self.super_root_mechanism = mechanism

    def set_remote_install(self, install):
        """Sets the type of Cantera installation

        Parameters
        ----------
        install : bool
            type of Cantera install

        """

        if install not in [True, False]:
            logger.warning("!!!WARNING!!! Argument not valid for remote_install, valid ones are :\
                   True (default) and False ")
            logger.warning('Default value will be kept')
        else:
            self.remote_install = install

    def set_overwrite(self, is_overwrite):
        """Sets if the cases should be overwritten or not

        Parameters
        ----------
        is_overwrite : bool
            if True, the data will be overwritten

        """

        self.overwrite = is_overwrite

    def set_remove_from(self, remove):
        """Sets if the invalid QSS species will be discarded from the list or the mechanism

        Parameters
        ----------
        remove : str
            'list' or 'mechanism'

        """

        if remove not in ['list', 'mechanism']:
            logger.warning("!!!WARNING!!! Argument not valid for remove_from, valid ones are :\
                   'list' and 'mechanism' (default)")
            logger.warning('Default value will be kept')
        else:
            self.remove_from = remove

    def set_diluent(self, diluent_list):
        """Sets a list of diluent species that will be kept even if they are not present in a reaction

        Parameters
        ----------
        diluent_list : list
            list of diluent species

        """

        if isinstance(diluent_list, str):
            if diluent_list not in self.mechanism.ctmech.species_names:
                logger.warning(diluent_list + ' not in mechanism, default N2 species will be set')
            else:
                self.diluent = [diluent_list]
        else:
            self.diluent = []
            for spec in diluent_list:
                if spec not in self.mechanism.ctmech.species_names:
                    logger.warning(spec + 'not in mechanism, this species has not been set as diluent')
                    self.diluent.append(spec)

    def set_stop_run(self, error_value, relative=False):
        """Sets variables for stop_value function

        Parameters
        ----------
        error_value : float
            threshold value of the error
        relative : bool
            if True, the threshold will be the input value times the maximum error allowed on the
            reduction (Default value = False)

        """

        self.error_threshold = error_value
        self.relative_threshold = relative

    def stop_value(self, max_error_tol):
        """Computes the value at which the run of the cases should stop

        Parameters
        ----------
        max_error_tol : flat
            maximum error tolerances for the reduction

        Returns
        -------
        stop_value : str
            value at which the run stops

        """

        if self.relative_threshold:
            stop_value = [self.error_threshold * err for err in max_error_tol]
        else:
            stop_value = [self.error_threshold for err in max_error_tol]

        return stop_value

    def set_qss_identification(self, method):
        """Switching QSS identification method

        Parameters
        ----------
        method : str
            if 'fast', reduction will use AOI method do identify QSS species

        """

        self.fast_qss_identification = bool(method in ['fast', 'AOI'])

    def set_restore(self, restore):
        """Set if cases are restored from existing data

        Parameters
        ----------
        restore : str
            set the restore method for `cases.run_cases()` method.
            Either 'auto', 'self' or 'parent'

        """

        self.restore = restore

    def set_reduction_outputs_dir(self, reduction_directory):
        """Output directory of the reduction data

        Parameters
        ----------
        reduction_directory : str
            reduction directory for outputs

        """

        self.reddir = reduction_directory
    
    def reduction_ranking(self, reduction_type, current_mechanism, sensitivity, has_qss, skeletal_mechanism):
        """Calculates the ranking coefficient for the reduction type given

        Parameters
        ----------
        reduction_type:
            reduction type, can either 'species', 'reactions', 'qss' 
        current_mechanism: 
            class `ARCANE.Mechanism` object
        sensitivity:
            specified sensitivity 
        has_qss:
            True if the mechanism contains already QSS species
        skeletal_mechanism:
            class `ARCANE.Mechanism` object corresponding to the skeletal mechanism if has_qss is True
        """
        # Apply DRGEP - Species to get ranked list of species
        if reduction_type in ['species', 'S']:
            species_ranked_dict = drgep.drgep_ranked(self.cases_list, current_mechanism, 'S', integrity=False,
                                                     sensitivity=sensitivity)
            # Initalize variables
            reduction_type = 'species'
            self.relevant_var = 'ns'
            self.name_of_keys = 'species'
            self.name_of_method = 'DRGEP'
            self.type_index = 0
            max_error_tol = self.max_error_spec
            self.ranked_dict = species_ranked_dict
            error_factor = self.error_factors[0]

        # Apply DRGEP - reactions to get ranked list of reactions
        elif reduction_type in ['reactions', 'R']:
            reactions_ranked_dict = drgep.drgep_ranked(self.cases_list, current_mechanism, 'R', sensitivity=sensitivity)

            # Initialize variables
            reduction_type = 'reactions'
            self.relevant_var = 'nr'
            self.name_of_keys = 'reactions'
            self.name_of_method = 'DRGEP'
            self.type_index = 1
            max_error_tol = self.max_error_reac
            self.ranked_dict = reactions_ranked_dict
            error_factor = self.error_factors[1]

        elif reduction_type in ['QSS', 'qss']:

            # Initalize variables
            reduction_type = 'QSS'
            self.relevant_var = 'nqss'
            self.name_of_keys = 'species'
            self.type_index = 2
            max_error_tol = self.max_error_qss
            error_factor = self.error_factors[3]

            if has_qss or not self.testing_reverse:
                testing_reverse = False
            else:
                testing_reverse = True

            while testing_reverse:

                # Testing the skeletal_mechanism against itself in f90 format
                fortran_skeletal_mechanism = mechanisms.Mechanism(skeletal_mechanism.path,
                                                                  f90='force_write', kinetics='custom',
                                                                  fit_order=self.fit_order)

                cases.run_cases(self.cases_list, fortran_skeletal_mechanism, overwrite=True)

                logger.info('\nError resulting from the fortran writing (super root vs skeletal with fortran kinetics)')
                error_values, error_super_list = error.error_dictionary(self.cases_list, self.super_root_mechanism,
                                                                        self.relevant_var, fortran_skeletal_mechanism,
                                                                        self.max_error_qss, error_factor=0.99)

                text_warning = """!!! WARNING !!! The error induced by the fortran is abnormally high.
                
The first reason for this error would be the presence in the mechanism of reactions with types other than:
reaction, three_body_reaction, falloff_reaction
If the pdep_reaction are present in your .cti, you can convert it with 
new_mechanism = old_mechanism.pdep_to_elementary()

If that's not the case, the error comes from the conversion of reversible reactions into irreversible reactions.
For this the backward rate constant is fitted over a range of temperature and if this fit is not good enough,           
An error can appear.

How to solve this problem ?

The quality of the fit is controlled by the discretisation of the temperatures range on which the constants are computed
The parameters are attributes of the Mechanism object: start_temperature, end_temperature, number_of_temperatures

A higher order fit (currently up to 4) or a function different from Arrhenius law could also be a solution.
For that case some coding could be required in order to solve the problem.

You should make a second script with only the comparison between classical and fortran kinetics,
And change the parameters (mainly end_temperature) until it reaches less than 1%
Here is an example:

fortran_skeletal_mechanism = mechanisms.Mechanism(skeletal_mechanism.path, f90='force_write', kinetics='custom',
                                                  start_temperature=300, end_temperature=3000, 
                                                  number_of_temperatures=6000)

cases.run_cases(cases_list, skeletal_mechanism)
cases.run_cases(cases_list, fortran_skeletal_mechanism, overwrite=True)

for case in cases_list:
    error.case_error(case, skeletal_mechanism, mechanism=fortran_skeletal_mechanism)

As this error is not frequent and requires special attention, if you still struggle, ask cazeres@cerfacs.fr 

"""

                if not all(error_super_list[index_err_v] < self.max_error_qss[index_err_v] * 0.99
                           for index_err_v, err_v in enumerate(error_super_list)):

                    text = f"!!! WARNING !!! The error induced by the fortran is abnormally high." \
                           f"A fit with a higher order Arrhenius will be done (order = {self.fit_order +1}) "

                    if self.fit_order < 4:
                        logger.warning(text)
                        self.fit_order += 1
                    else:
                        logger.warning(text_warning)
                        quit()

                else:
                    testing_reverse = False

            # If root mechanism has QSS species, the LOI method is not useable
            if not self.fast_qss_identification and not has_qss:
                # Apply LOI to get ranked list of QSS species
                qss_ranked_dict = qss.qss_ranked_species_LOI(self.cases_list, current_mechanism)
                self.name_of_method = 'LOI'
            else:
                # Apply AOI to get ranked list of QSS species
                qss_ranked_dict = qss.qss_ranked_species_AOI(self.cases_list, current_mechanism)
                self.name_of_method = 'AOI'

            self.ranked_dict = qss_ranked_dict

        else:
            self.relevant_var = ''
            max_error_tol = 0
            self.ranked_dict = {}
            logger.error('This reduction type does not exist !')
            quit()
        
        return reduction_type, max_error_tol, error_factor

    def initialize_max_errors(self):
        """Initialize the maximum errors list

        """
        for case in self.cases_list:
            for err_method in case.error_dict:
                # Store the error method for the species variable
                if err_method in kwdict.names['Y']:
                    for target in case.targets:
                        if target != 'HeatRelease':
                            self.max_errors.append(case.error_dict[err_method])
                            self.max_error_spec.append(case.error_dict[err_method] * self.error_factors[0])
                            self.max_error_reac.append(case.error_dict[err_method] * self.error_factors[1])
                            self.max_error_qss.append(case.error_dict[err_method] * self.error_factors[3])

                # Store the error for the global variable
                else:
                    self.max_errors.append(case.error_dict[err_method])
                    self.max_error_spec.append(case.error_dict[err_method] * self.error_factors[0])
                    self.max_error_reac.append(case.error_dict[err_method] * self.error_factors[1])
                    self.max_error_qss.append(case.error_dict[err_method] * self.error_factors[3])
        
    def reduction(self, reduction_type, sensitivity=0, permissive=False):
        """
        Performs the reduction of a mechanism base on a list of cases

        Parameters
        ----------
        reduction_type : str
            type of reduction ('species' or 'reactions')
        sensitivity : bool
            activates the DRGEP-ASA --> will be available in v1.1 (Default value = 0)
        permissive : bool
            says if case allows failing cases or not (Default value = False)

        Returns
        -------
        reduced_mechanism :
            reduced :func:`~ARCANE.mechanisms.Mechanism` object

        """
        ##################
        # Initialization #
        ##################

        # Initialize header
        head, tail = get_head_and_tail(reduction_type)

        # Shortcut allowing to use reduction for the lumping reduction
        if reduction_type in ['lumping', 'L']:
            reduced_mechanism = self.lumping_reduction()
            return reduced_mechanism

        # Initialize case directory
        database.create_dir(self.reddir, False)

        # Useful mechanism instantiation
        root_mechanism = self.mechanism
        super_root_mechanism = self.super_root_mechanism

        logger.goodnews("\nThe super root mechanism is " + super_root_mechanism.nickname)
        logger.goodnews("\nThe root mechanism is " + root_mechanism.nickname + "\n")

        # Finds if root mechanism has QSS species
        has_qss, skeletal_mechanism, species_qss_names =  root_mechanism.has_qss()
        species_qss_names_all = species_qss_names[:]
        species_qss_names_comp = species_qss_names[:]
 
        current_mechanism = root_mechanism
        reduced_mechanism = current_mechanism

        # Setting diluent species
        root_mechanism.inert = self.diluent
        current_mechanism.inert = self.diluent
        reduced_mechanism.inert = self.diluent

        # Since the handling of pdep reactions in the f90 is not good,
        # Converting all to elementary reactions with the 1 bar arrhenius
        if reduction_type in ['QSS', 'qss'] and self.pdep_value:
            skeletal_mechanism = skeletal_mechanism.pdep_to_elementary(pressure=self.pdep_value)

        # If a super root is included, include it in the list
        if super_root_mechanism and super_root_mechanism != root_mechanism:
            mechanisms_list = [super_root_mechanism, current_mechanism]
            root_index = 1
        else:
            mechanisms_list = [current_mechanism]
            root_index = 0

        # Initialize database
        mechanism_properties = pd.DataFrame(data=[], index=['previous', 'species', 'reactions',
                                                            'qss', 'Error maximum', 'Error maximum variable', 'fail'])
        mechanism_properties[root_mechanism.name] = ""
        mechanism_properties[root_mechanism.name]['previous'] = root_mechanism.name
        mechanism_properties[root_mechanism.name]['species'] = root_mechanism.ctmech.species_names
        mechanism_properties[root_mechanism.name]['reactions'] = root_mechanism.ctmech.reaction_equations()

        mechanism_properties[root_mechanism.name]['Error maximum'] = 0
        mechanism_properties[root_mechanism.name]['Error maximum variable'] = 'tig'

        # Recover cases_list
        cases_list = self.cases_list

        # Validation mechanism flags
        valid_mechanism = False
        one_valid_mechanism = False
        found_reduced = False
        fails = 0

        # Error directories
        store_error = {}
        store_dictionary = {}

        # Maximum error tolerance
        self.max_errors = []
        self.max_error_spec = []
        self.max_error_reac = []
        self.max_error_qss = []

        # Flag for first passage
        first_passage = 0
        first_passage_non_stop = 0

        # Setting non-stop parameter
        if isinstance(self.non_stop, bool):
            self.non_stop = [self.non_stop] * 3

        # Setting dichotomy parameters
        nb_big = current_mechanism.ns
        nb_small = 0

        # Setting default values
        if self.error_factors == [0, 0, 0, 0]:
            self.error_factors = [1, 1, 1, 1]

        # Initializing error vectors
        self.initialize_max_errors()
 
        ##############################
        # Run cases before reduction #
        ##############################

        # Retrieving important species that cannot be discarded
        all_targets = []
        for case in cases_list:
            for target in case.targets:
                if target != 'HeatRelease':
                    all_targets.append(target)

        # Run cases
        self.cantera_run_time += cases.run_cases(cases_list, super_root_mechanism, permissive=permissive, timing=True)
        self.cantera_run_time += cases.run_cases(cases_list, current_mechanism, restore=self.restore,
                                                 permissive=permissive, timing=True)

        # Resets the dynamic library
        current_mechanism.reset()

        # Stops here if the reference cases are crashing
        one_crash = False
        for case in cases_list:
            if not case.success:
                logger.error(case.myid + ' simulation failed')
                one_crash = True
        if one_crash:
            logger.error('ERROR ! Process stopped there because of references cases failure')
            quit()

        #############################
        # Relevant variable ranking #
        #############################

        # Ranking functions
        reduction_type, max_error_tol, error_factor = self.reduction_ranking(reduction_type, current_mechanism, 
                                                                             sensitivity, has_qss, skeletal_mechanism)

        #####################
        # Step of reduction #
        #####################

        if self.step == 'auto':
            step_value = 1
            if self.relevant_var in ['ns', 'nqss']:
                for threshold, step in self.ns_auto_stepping:
                    if current_mechanism.ns > threshold:
                        step_value = step
                        break
            else:
                for threshold, step in self.nr_auto_stepping:
                    if current_mechanism.nr > threshold:
                        step_value = step
                        break

        elif self.dichotom:
            if self.relevant_var in ['ns', 'nqss']:
                step_value = int(root_mechanism.ns / 2)
            else:
                step_value = int(root_mechanism.nr / 2)
        else:
            step_value = int(self.step)

        # Further manipulation of ranking

        # Setting diluent species to 1
        for key in self.ranked_dict:
            if key in self.diluent:
                self.ranked_dict[key] = 1
        # Clipping values that are too low
            if self.ranked_dict[key] < 1e-60:
                self.ranked_dict[key] = 0.0

        # Setting the value of the ranked dictionary as too high if stop_criterion is activated
        if self.stop_criterion and self.non_stop[self.type_index]:
            too_high_value = tools.cutting_criterion(self.ranked_dict) * self.criterion_factors[self.type_index]
            if too_high_value > 1:
                too_high_value = 1
        else:
            if sensitivity:
                too_high_value = 10 ** 50
            else:
                too_high_value = 1

        # Displaying the list of sorted items with their coefficient
        logger.info("Sorted list of " + self.name_of_keys + " and their " + self.name_of_method + " coefficients\n")
        sorted_list = sorted(zip(self.ranked_dict.values(), self.ranked_dict.keys()), key=lambda x: x[0])
        line_plotted = False
        for s in sorted_list:
            # Drawing a line separating the too high values
            if s[0] > too_high_value and not line_plotted and not self.bypass_ones and not valid_mechanism:
                logger.info('-' * len(max([str(s[1]) for s in sorted_list], key=len)) + '----------------')
                line_plotted = True
            # Displaying relevant variable and coefficient
            logger.info(str(s[1]).ljust(len(max([str(s[1]) for s in sorted_list], key=len)))
                        + " => " + str('{:.6E}'.format(s[0])))

        # Initialization of plotting data
        ranked = sorted(zip(self.ranked_dict.keys(), self.ranked_dict.values()), key=lambda x: x[1])

        # Automatically discard species/reactions with null coefficients
        to_discard = []
        start = 0
        if self.skip_zeros:
            while ranked[start][1] == 0:
                to_discard.append(ranked[start][0])
                start = start + 1
        index = start

        #######################
        # Error instantiation #
        #######################

        # Stop value fixed to the error if non_stop
        if self.non_stop[self.type_index] or reduction_type == 'R':
            self.set_stop_run(1, relative=True)
        stop_value = self.stop_value(max_error_tol)

        store_error[getattr(current_mechanism, self.relevant_var)] = np.zeros([len(max_error_tol)])

        ##############################
        # Reduction of the mechanism #
        ##############################
        while index < len(ranked):

            # Stop directly if the rank is over the cutting criterion or if the
            # ranked items equal to 1 are not bypassed
            if ranked[index][1] >= too_high_value and not self.bypass_ones \
                    or ranked[index][0] in all_targets:
                logger.info('The following values of the ranked items are considered too high for the reduction to work'
                            '\n if you still want to compute them, set the stop_criterion to False\n')
                if not found_reduced:
                    reduced_mechanism = mechanisms_list[-1]
                break
            else:
                # parent_mechanism is the father of the reduced mechanism
                parent_mechanism = current_mechanism

                ###################################################
                # Items to discard and reduced mechanism creation #
                ###################################################

                # Restart list if not first passage.
                # Otherwise, keep preceding list.
                if first_passage != 0:
                    to_discard = []

                # Case reduction type is species one
                if reduction_type in ['species', 'S']:
                    to_discard.extend([ranked[index + increment][0] for increment in range(step_value)
                                       if ranked[index + increment][0] not in self.diluent
                                       and ranked[index + increment][0] in
                                       current_mechanism.ctmech.species_names])
                    # If function is empty (because species cannot be taken for
                    # instance), go to next step.
                    if not to_discard:
                        index += step_value
                        continue

                    logger.info('\nRemoving species : ' + ', '.join(to_discard))

                    # Build reduced mechanism
                    current_mechanism = parent_mechanism.remove_species(to_discard)
                    current_mechanism.how = 'rmS'

                # Case reduction type is reactions one
                elif reduction_type in ['reactions', 'R']:
                    to_discard.extend([ranked[index + increment][0] for increment in range(step_value)
                                       if ranked[index + increment][0] not in
                                       self.diluent])

                    logger.info('\nRemoving reactions : ' + ' \n'.join([str(equation) for equation in to_discard]))

                    # Build reduced mechanism
                    current_mechanism = parent_mechanism.remove_reactions(to_discard)
                    current_mechanism.how = 'rmR'

                # Case reduction type is QSS species one
                elif reduction_type in ['qss', 'QSS']:
                    to_discard.extend([ranked[index + increment][0] for increment in range(step_value)
                                       if ranked[index + increment][0] not in
                                       self.diluent])
                    if first_passage > 0 and self.non_stop[self.type_index]:
                        if not valid_mechanism:
                            if len(species_qss_names_all) != 0:
                                species_qss_names_all = species_qss_names_all[:-1]
                            else:
                                species_qss_names_all = []

                    [species_qss_names_all.append(spec) for spec in to_discard]

                    # Build cantera solution
                    new_ctmech, species_qss_names = custom.problematic_qss(skeletal_mechanism,
                                                                           species_qss_names_all,
                                                                           remove_from=self.remove_from)

                    # Build skeletal mechanism
                    new_skeletal_mechanism = mechanisms.Mechanism(new_ctmech,
                                                                  parent=skeletal_mechanism,
                                                                  how='errorcheck',
                                                                  species_qss_names=[])

                    # Check to ensure that when using "remove_from = list", the
                    # number of QSS species is always increasing
                    if species_qss_names == species_qss_names_comp:
                        continue
                    elif len(species_qss_names) < len(species_qss_names_comp):
                        continue
                    else:
                        species_qss_names_comp = species_qss_names[:]

                    logger.info('\nNew set of QSS species : ' + ', '.join(species_qss_names))
                    if len(to_discard) < len(species_qss_names):
                        logger.info('(Species added to previous set : ' + ', '.join(to_discard) + ')')

                    # If overwrite, the f90 will also be rewritten
                    if self.overwrite:
                        write_rule = 'force_write'
                    else:
                        write_rule = 'write'

                    # Create Cantera solution for qss run (cti printed)
                    reduced_how = 'QSS'

                    current_mechanism = mechanisms.Mechanism(new_skeletal_mechanism.ctmech,
                                                             parent=parent_mechanism,
                                                             how=reduced_how,
                                                             species_qss_names=species_qss_names,
                                                             f90=write_rule,
                                                             fit_order=self.fit_order,
                                                             start_temperature=self.start_temperature,
                                                             end_temperature=self.end_temperature,
                                                             number_of_temperatures=self.number_of_temperatures,
                                                             fit_pressure=self.fit_pressure)
                else:
                    logger.error('The reduction type ' + reduction_type + ' does not exist !')
                    quit()

                mechanisms_list.append(current_mechanism)

                ###############################
                # Run reduced mechanism cases #
                ###############################

                self.cantera_run_time += cases.run_cases(cases_list, current_mechanism, restore=self.restore,
                                                         overwrite=self.overwrite, permissive=permissive, timing=True)

                one_case_failed = False
                for case in cases_list:
                    if not case.success:
                        one_case_failed = True

                if one_case_failed and not self.non_stop[self.type_index]:
                    logger.warning('WARNING: One case failed so this reduction step will stop there')
                    break

                # Resets dynamic library
                current_mechanism.reset()

                ####################################
                # Calculate the error and store it #
                ####################################

                # Error determination
                if super_root_mechanism:
                    error_values, error_super_list = error.error_dictionary(self.cases_list, super_root_mechanism,
                                                                            self.relevant_var, current_mechanism,
                                                                            max_error_tol, error_factor=error_factor)
                else:
                    error_values, error_super_list = error.error_dictionary(self.cases_list, root_mechanism,
                                                                            self.relevant_var, current_mechanism,
                                                                            max_error_tol, error_factor=error_factor)

                # Store important information
                mechanism_properties[current_mechanism.name] = ""
                mechanism_properties[current_mechanism.name]['previous'] = parent_mechanism.name
                mechanism_properties[current_mechanism.name]['species'] = current_mechanism.ctmech.species_names
                mechanism_properties[current_mechanism.name]['reactions'] \
                    = current_mechanism.ctmech.reaction_equations()
                mechanism_properties[current_mechanism.name]['qss'] = species_qss_names

                mechanism_properties[current_mechanism.name]['Error maximum'] = np.max(error_super_list)
                mechanism_properties[current_mechanism.name]['Error maximum variable'] = 'tig'

                ###############################################
                # Check to tell if the mechanism is OK or not #
                ###############################################

                # Mechanism is not valid at the beginning
                valid_mechanism = False

                # Compare errors to the errors defined by the user and state on
                # the validity of the mechanism.
                if all(error_super_list[index_err_v] < max_error_tol[index_err_v]
                       for index_err_v, err_v in enumerate(error_super_list)):
                    valid_mechanism = True
                    one_valid_mechanism = True
                    mechanism_properties[mechanisms_list[-1].name]['fail'] = 'no'
                    logger.debug('The mechanism ' + mechanisms_list[-1].name + ' is valid.')

                    # Removing previous mechanisms data from memory
                    if parent_mechanism.name not in [root_mechanism.name, super_root_mechanism.name]:
                        for case in self.cases_list:
                            del case.stored_data[parent_mechanism.name]

                    # Number of fails comes back to zero
                    if fails >= 1:
                        fails = 0
                else:
                    # Number of fails in increased
                    fails += 1
                    mechanism_properties[mechanisms_list[-1].name]['fail'] = 'max'

                    # Removing current mechanism data from memory
                    if parent_mechanism.name not in [root_mechanism.name, super_root_mechanism.name]:
                        for case in self.cases_list:
                            if current_mechanism.name in case.stored_data:
                                del case.stored_data[current_mechanism.name]

                # Flag on whether the first mechanism where there is an error
                # is chosen or the last one
                if self.which_reduced == 'first':
                    # If the first mechanism is chosen, when there is an error,
                    # the mechanism found is the preceding one.
                    if not all(error_super_list[index_err_v] < max_error_tol[index_err_v]
                               for index_err_v, err_v in enumerate(error_super_list)) \
                            and not found_reduced:
                        reduced_mechanism = mechanisms_list[-2]
                        found_reduced = True
                else:
                    # Otherwise, if there is compensation and that the maximum
                    # compensation number is reached, the mechanism come back
                    # to the compensation number before and fails return to
                    # zero
                    if fails == self.compensation and self.compensation != 0:
                        current_mechanism = mechanisms_list[-self.compensation - 1]
                        fails = 0
                        mechanisms_list = mechanisms_list[:-self.compensation - 1]
                    # In all cases, if the mechanism is valid, it is chosen and
                    # found as the reduced mechanism
                    if valid_mechanism:
                        reduced_mechanism = mechanisms_list[-1]
                        found_reduced = True

                # Check if the error goes above one of the stop_value
                if not (all(error_super_list[index_err_v] <= stop_value[index_err_v]
                            for index_err_v, err_v in enumerate(error_super_list))):
                    mechanism_properties[mechanisms_list[-1].name]['fail'] = 'stop'

                # Store the error for case considered (species, equation or
                # qss) as well as data for further plots except when non stop
                # fails
                store_error[current_mechanism.name] = []
                if self.non_stop[self.type_index] and mechanism_properties[current_mechanism.name]['fail'] == 'stop':
                    pass
                else:
                    # Store data for plotting
                    max_error = []
                    av_error = []
                    max_tol = []
                    l_tol = 0

                    for x in error_values.values():
                        for y in x:
                            store_error[current_mechanism.name].append(y)

                    store_dictionary[current_mechanism.name] = error_values

                    # Store maximum and average plot containing the same type of case, quantity, method and tolerance
                    for idx, x in enumerate(store_dictionary[current_mechanism.name].values()):
                        max_error.append(np.max(x))
                        av_error.append(np.mean(x))
                        max_tol.append(max_error_tol[l_tol + len(x) - 1])
                        l_tol = len(x) + l_tol

                    # Initiate vector at first passage to plot the maximum and average of the error
                    # for all the type of errors
                    if first_passage_non_stop == 0:
                        self.xdata = []
                        ydata_max = [[] for x in range(len(max_error))]
                        ydata_av = [[] for x in range(len(max_error))]
                        first_passage_non_stop = 1

                    # Create the vectors of error that need to be plotted
                    self.xdata.append(getattr(current_mechanism, self.relevant_var))

                    for index_y, y in enumerate(max_error):
                        ydata_max[index_y].append(y)
                    for index_y, y in enumerate(av_error):
                        ydata_av[index_y].append(y)

                # Creation of the directed graph
                # (must be done before exiting the code because error might be
                # bigger than stop value)
                dot_text = """
digraph {
    graph [splines=ortho]
    node [shape=record];
"""
                n_split = self.reduction_graph
                if n_split > 0:
                    if reduction_type in ['species', 'S']:
                        argument = 'species'
                    elif reduction_type in ['reactions', 'R']:
                        argument = 'reactions'
                    elif reduction_type in ['QSS', 'qss']:
                        argument = 'qss'
                    else:
                        logger.error('Reduction type is not good !')

                    dot_text += "{ rank=same; " + root_mechanism.name + " "
                    for n in range(1, mechanism_properties.shape[1], 1):
                        if len(mechanism_properties[mechanism_properties.keys()[n]][argument]) % n_split == 0:
                            dot_text += mechanism_properties.keys()[n] + " "
                    dot_text += " }\n"

                    for n in range(1, mechanism_properties.shape[1], 1):
                        dot_text += mechanism_properties[mechanism_properties.keys()[n]]['previous'] + \
                                    " -> " + mechanism_properties.keys()[n]

                        dot_text += ' [label="'
                        previous_n = mechanism_properties[mechanism_properties.keys()[n]]['previous']
                        set_n = mechanism_properties[previous_n][argument]
                        set_n1 = mechanism_properties[mechanism_properties.keys()[n]][argument]
                        for element in set(set_n).symmetric_difference(set(set_n1)):
                            dot_text += " - " + str(element)
                        dot_text += '",'
                        if mechanism_properties[mechanism_properties.keys()[n]]['fail'] == 'no':
                            dot_text += 'color=green,'
                        elif mechanism_properties[mechanism_properties.keys()[n]]['fail'] == 'max':
                            dot_text += 'color=orange,'
                        else:
                            dot_text += 'color=red,'
                        if len(mechanism_properties[mechanism_properties.keys()[n]][argument]) % n_split == 0:
                            dot_text += 'weight="1"];\n'
                        else:
                            dot_text += 'weight="5"];\n'
                    dot_text += "}"

                # Says if it is first passage or not
                if first_passage == 0:
                    first_passage = 1

                # Instructions to do in case of stopping
                # - In case of non stop and no compensation, the code is
                # continuing with the preceding mechanism.
                # - In case of compensation, nothing appends.
                # - In all other cases, the code breaks the loop.
                if mechanism_properties[mechanisms_list[-1].name]['fail'] == 'stop':
                    if self.non_stop[self.type_index] and self.compensation == 0:
                        current_mechanism = reduced_mechanism
                    elif self.compensation or self.dichotom:
                        pass
                    else:
                        break

                ######################
                # Index modification #
                ######################
                if self.compensation != 0 and fails == self.compensation:
                    # If the value reaches the compensation number, the code comes
                    # back to the compensation minus 1 index
                    index = index - (self.compensation - 1)
                elif self.dichotom:
                    if valid_mechanism:
                        index = index + step_value
                else:
                    # In all other cases, it just increases from the step_value
                    index = index + step_value

                ######################
                # Step instantiation #
                ######################

                if self.step == 'auto':
                    step_value = 1
                    if self.relevant_var in ['ns', 'nqss']:
                        for threshold, step in self.ns_auto_stepping:
                            if current_mechanism.ns > threshold:
                                step_value = step
                                break
                    else:
                        for threshold, step in self.nr_auto_stepping:
                            if current_mechanism.nr > threshold:
                                step_value = step
                                break

                elif self.dichotom:
                    """
                    if valid_mechanism:
                        if self.relevant_var in ['ns', 'nqss']:
                            #nb_species_a = current_mechanism.ns
                            step_value = int(current_mechanism.ns / 2)
                        else:
                            step_value = int(current_mechanism.nr / 2)
                    else:
                        dicho_ns = current_mechanism.ns
                        dicho_nr = current_mechanism.nr
                        mechanisms_list.remove(current_mechanism)
                        current_mechanism = mechanisms_list[-1]
                        if self.relevant_var in ['ns', 'nqss']:
                            step_value = int(current_mechanism.ns - current_mechanism.ns / 2 - dicho_ns / 2)
                        else:
                            step_value = int(current_mechanism.nr - current_mechanism.nr / 2 - dicho_nr / 2)
                    """

                    if valid_mechanism:
                        if reduction_type in ['species', 'S', 'QSS', 'qss']:
                            nb_big = current_mechanism.ns
                        else:
                            nb_big = current_mechanism.nr
                    else:
                        if reduction_type in ['species', 'S', 'QSS', 'qss']:
                            nb_small = current_mechanism.ns
                        else:
                            nb_small = current_mechanism.nr

                    nb_mean = int((nb_big + nb_small) / 2)
                    if not valid_mechanism:
                        mechanisms_list.remove(current_mechanism)
                        current_mechanism = mechanisms_list[-1]

                    if reduction_type in ['species', 'S', 'QSS', 'qss']:
                        step_value = current_mechanism.ns - nb_mean
                    else:
                        step_value = current_mechanism.nr - nb_mean

                    if nb_big - nb_small <= 1:
                        break

        #####################
        # Dot graph writing #
        #####################

        if reduction_type not in ['lumping', 'L'] and self.reduction_graph != 0:
            i = 0
            while os.path.exists("graphs/dotgraph%s" % i):
                i += 1
            write_dot = open("graphs/dotgraph%s" % i, "w+")
            write_dot.write(dot_text)
            write_dot.close()
            dot_graph = Source.from_file('graphs/dotgraph%s' % i)
            dot_graph.render("graphs/dotgraph%s" % i, view=False, format='png')

        #####################################################################
        # Creation of the reduced mechanism if not created during reduction #
        #####################################################################

        # Might happen in the following cases :

        # Not one single mechanism was OK at the end of the reduction. Then the
        # reduced mechanism is the root mechanism.
        if not one_valid_mechanism:
            reduced_mechanism = mechanisms_list[root_index]

        reduced_mechanism.copy(new_name='reduced' + reduced_mechanism.name, directory=self.reddir)

        ########################
        # Print useful outputs #
        ########################

        # Message telling the reduction is finished
        head, tail = display.head_tail_charac(text='Reduction finished', size='small', style='@')
        logger.goodnews(head)

        # Building the graph error = f(relevant_var)
        if reduced_mechanism.name != root_mechanism.name:
            self.assess_reduction(reduction_type, reduced_mechanism, root_mechanism,
                                  store_dictionary, error_values, ydata_max, ydata_av, max_tol, max_error_tol)

        else:
            logger.goodnews('Optimal mechanism is the root mechanism ' + str(reduced_mechanism.name))
            logger.goodnews(tail)

        # Store in database with different name for reduced mechanism
        reduced_mechanism.copy(new_name=f'Red_{reduced_mechanism.name}')

        return reduced_mechanism

    def lumping_reduction(self, lumping_steps='all', permissive=False):
        """
        Performs the reduction of a mechanism base on a list of cases

        Parameters
        ----------
        lumping_steps : str
            name of the lumping step ('all', 'identical', 'brute', '2') (Default value = 'all')
        permissive : bool
            boolean to allow for further reduction (Default value = False)

        Returns
        -------
        reduced_mechanism :
            reduced :func:`~ARCANE.mechanisms.Mechanism` object

        """

        head, tail = display.head_tail_charac(text='Lumping reduction', size='big', style='#')

        logger.info(head)

        # Initialize case directory
        database.create_dir(self.reddir, False)

        root_mechanism = self.mechanism
        super_root_mechanism = self.super_root_mechanism
        current_mechanism = root_mechanism
        reduced_mechanism = current_mechanism

        cases_list = self.cases_list

        one_valid_mechanism = False

        # Maximum error tolerance
        max_error_lumping = []
        self.max_errors = []

        # Setting default values
        if self.error_factors == [0, 0, 0, 0]:
            self.error_factors = [1, 1, 1, 1]

        # Retrieving important species that cannot be discarded
        all_targets = []
        for case in cases_list:
            if hasattr(case, 'targets'):
                for target in case.targets:
                    if target != 'HeatRelease':
                        all_targets.append(target)

        for case in cases_list:
            for err_method in case.error_dict:
                # Store the error method for the species variable
                if err_method in kwdict.names['Y']:
                    for target in case.targets:
                        if target != 'HeatRelease':
                            self.max_errors.append(case.error_dict[err_method])
                            max_error_lumping.append(case.error_dict[err_method] * self.error_factors[2])

                # Store the error for the global variable
                else:
                    self.max_errors.append(case.error_dict[err_method])
                    max_error_lumping.append(case.error_dict[err_method] * self.error_factors[2])

        # Run cases
        self.cantera_run_time += cases.run_cases(cases_list, super_root_mechanism, permissive=permissive, timing=True)
        self.cantera_run_time += cases.run_cases(cases_list, current_mechanism, restore=self.restore,
                                                 overwrite=self.overwrite, permissive=permissive, timing=True)

        # Stops here if the reference cases are crashing
        one_crash = False
        for case in cases_list:
            if not case.success:
                logger.error(case.myid + ' simulation failed')
                one_crash = True
        if one_crash:
            logger.error('ERROR ! Process stopped there because of references cases failure')
            quit()

        # Error directories
        store_error = {}
        store_dictionary = {}

        if lumping_steps == 'all':
            lumping_steps = ['Identical species lumping', 'Brute force lumping', '2 by 2 lumping']
        elif lumping_steps == 'identical':
            lumping_steps = ['Identical species lumping']
        elif lumping_steps == 'brute':
            lumping_steps = ['Brute force lumping']
        elif lumping_steps == '2':
            lumping_steps = ['2 by 2 lumping']
        else:
            logger.error('ERROR ! Lumping step name is not valid')

        self.relevant_var = 'ns'
        max_error_tol = max_error_lumping

        first_passage = 0

        # stop value fixed to the error if non_stop
        self.set_stop_run(1, relative=True)
        stop_value = self.stop_value(max_error_tol)

        store_error[getattr(current_mechanism, self.relevant_var)] = np.zeros([len(max_error_tol)])

        # Set it as current mechanism (no duplicates here, same object)
        current_mechanism = root_mechanism

        # If a super root is included, include it in the list
        if super_root_mechanism and super_root_mechanism != root_mechanism:
            mechanisms_list = [super_root_mechanism, current_mechanism]
            root_index = 1
        else:
            mechanisms_list = [current_mechanism]
            root_index = 0

        for lumping_step in lumping_steps:

            head, tail = display.head_tail_charac(text=lumping_step, size='medium', style='#')

            logger.info(head)

            isomers_dict = lumping.find_isomers(current_mechanism)

            species_to_lump_list = list(isomers_dict.values())
            species_to_lump_list.reverse()

            if lumping_step == 'Identical species lumping':

                identical_species = []
                for isomers in species_to_lump_list:
                    identical_species += lumping.detect_identical_species(current_mechanism, isomers)

                ranked = identical_species

            elif lumping_step == 'Brute force lumping':

                ranked = species_to_lump_list

            else:
                lumpable_species = []
                for isomers in species_to_lump_list:
                    temp_list, species_ratio = lumping.find_species_ratios(cases_list, current_mechanism, isomers)
                    if temp_list:
                        spec_arrhenius_dict = lumping.create_spec_arrhenius_dict(temp_list, species_ratio)
                    else:
                        break

                    rough_groups = lumping.lump_species_2_by_2(current_mechanism, spec_arrhenius_dict, 0.5)

                    if rough_groups:
                        lumpable_species += lumping.get_best_combination(rough_groups)[0]

                ranked = lumpable_species

            if not ranked:
                logger.info('No lumpable species were found in this step\n')

            for index, to_lump_list in enumerate(ranked):

                parent_mechanism = current_mechanism

                # Unable to lump inlet species at the moment
                inlet_species = all_targets
                for case in cases_list:
                    if hasattr(case, 'composition'):
                        inlet_species += list(case.composition.keys())

                to_lump_list = [spec for spec in to_lump_list if spec not in inlet_species]
                if len(to_lump_list) <= 1:
                    continue
                else:
                    temp_list, species_ratio = lumping.find_species_ratios(cases_list, current_mechanism, to_lump_list)
                    if temp_list:
                        spec_arrhenius_dict = lumping.create_spec_arrhenius_dict(temp_list, species_ratio)
                        logger.info('\nLumping species : ' + ', '.join(to_lump_list))

                        reduced_ctmech = lumping.lump_species(current_mechanism, spec_arrhenius_dict)
                    else:
                        break

                    # Append new mechanism
                    current_mechanism = mechanisms.Mechanism(parent=parent_mechanism,
                                                             cti=reduced_ctmech,
                                                             how='lumping')

                    mechanisms_list.append(current_mechanism)

                    # Run cases
                    self.cantera_run_time += cases.run_cases(cases_list, current_mechanism, restore=self.restore,
                                                             overwrite=self.overwrite, permissive=permissive,
                                                             timing=True)

                    # Error determination
                    if super_root_mechanism:
                        error_values, error_super_list = error.error_dictionary(self.cases_list, super_root_mechanism,
                                                                                self.relevant_var, current_mechanism,
                                                                                max_error_tol)
                    else:
                        error_values, error_super_list = error.error_dictionary(self.cases_list, root_mechanism,
                                                                                self.relevant_var, current_mechanism,
                                                                                max_error_tol)

                    # Store the error for case considered (species, equation or qss)
                    store_error[current_mechanism.name] = []
                    for x in error_values.values():
                        for y in x:
                            store_error[current_mechanism.name].append(y)

                    store_dictionary[current_mechanism.name] = error_values

                    # Store data for plotting
                    max_error = []
                    av_error = []
                    max_tol = []
                    l_tol = 0

                    # Store maximum and average plot containing the same type of case, quantity, method and tolerance
                    for idx, x in enumerate(error_values.values()):
                        max_error.append(np.max(x))
                        av_error.append(np.mean(x))
                        max_tol.append(max_error_tol[l_tol + len(x) - 1])
                        l_tol = len(x) + l_tol

                    # Initiate vector at first passage to plot the maximum and average of the error
                    # for all the type of errors
                    if first_passage == 0:
                        self.xdata = []
                        ydata_max = [[] for x in range(len(max_error))]
                        ydata_av = [[] for x in range(len(max_error))]
                        first_passage = 1

                    # Create the vectors of error that need to be plotted
                    self.xdata.append(getattr(current_mechanism, self.relevant_var))

                    for index_y, y in enumerate(max_error):
                        ydata_max[index_y].append(y)
                    for index_y, y in enumerate(av_error):
                        ydata_av[index_y].append(y)

                    if all(error_super_list[index_err_v] < max_error_tol[index_err_v]
                           for index_err_v, err_v in enumerate(error_super_list)):
                        reduced_mechanism = mechanisms_list[-1]
                        one_valid_mechanism = True

                    current_mechanism = reduced_mechanism

            if not one_valid_mechanism:
                reduced_mechanism = mechanisms_list[root_index]

            reduced_mechanism.copy(new_name='reduced' + reduced_mechanism.name, directory=self.reddir)

            # Print useful outputs and graphs error = f(relevant_var)
            if lumping_step == lumping_steps[-1]:
                head, tail = display.head_tail_charac(text='Reduction finished', size='small', style='@')
            else:
                head, tail = display.head_tail_charac(text='Step finished', size='small', style='@')
            logger.goodnews(head)

            if reduced_mechanism.name != root_mechanism.name:
                self.assess_reduction('lumping', reduced_mechanism, root_mechanism, store_dictionary,
                                      error_values, ydata_max, ydata_av, max_tol, max_error_tol)
            else:
                logger.goodnews('Optimal mechanism is the root mechanism ' + str(reduced_mechanism.name))
                logger.goodnews(tail)

        # Store in database with different name for reduced mechanism
        reduced_mechanism.copy(new_name=f'Red_{reduced_mechanism.name}')

        return reduced_mechanism

    def assess_reduction(self, reduction_type, reduced_mechanism, root_mechanism,
                         store_dictionary, error_values, ydata_max, ydata_av, max_tol, max_error_tol):
        """
        Build the error graphs and calculate the maximum and average errors

        Parameters
        ----------
        reduction_type : str
            string to specify the type of reduction
        reduced_mechanism : str
            string to specify the reduced mechanism
        root_mechanism : str
            string to specify the root mechanism
        relevant_var : str
            relevant variable associated to the reduction_type
        store_dictionary : dict
            dictionary where the errors are stored by cases
        error_values : dict
            dictionary where the errors are stored linearly
        ydata_max : list
            vector with the maximum errors of values
        ydata_av : list
            vector with the average errors of values
        max_tol : list
            vector of the maximum tolerances for the errors
        max_error_tol : float
            maximum of the error tolerances

        """

        head, tail = display.head_tail_charac(text='', size='small', style='@')

        # Store the maximum and average errors of the mechanism
        max_error = []
        av_error = []
        for x in store_dictionary[reduced_mechanism.name].values():
            max_error.append(np.max(x))
            av_error.append(np.mean(x))

        # Print interesting data
        logger.goodnews('Optimal mechanism with ' + str(self.relevant_var) + " = " + str(getattr(reduced_mechanism,
                                                                                            self.relevant_var)))
        for index_x, x in enumerate(store_dictionary[reduced_mechanism.name]):

            reactor_type_name = kwdict.reactor_labels[x[1]][0]

            error_name = kwdict.get_full_name(x[2])

            if callable(x[2]):
                error_name = x[2].__name__
            method_name = x[3]
            tolerence_string = x[4]

            separator = ', '

            if method_name == 'none':
                method_name = ''
                separator = ''

            logger.info('Maximum error on ' + reactor_type_name + ', ' + error_name + separator
                        + method_name + ', ' + tolerence_string + ' %  = ' +
                        '{:.2f}'.format(max_error[index_x] * 100) + ' %')
            logger.info('Average error on ' + reactor_type_name + ', ' + error_name + separator
                        + method_name + ', ' + tolerence_string + ' %  = ' +
                        '{:.2f}'.format(av_error[index_x] * 100) + ' %')
        logger.goodnews('!!! INFO !!! .cti file for the reduced mechanism has been written as ' +
                    'reduced' + reduced_mechanism.name + '.cti')

        logger.goodnews(tail)

        # Create graphs and register them
        if self.error_plot and len(self.xdata) >= 2 and getattr(root_mechanism, self.relevant_var) != \
                getattr(reduced_mechanism, self.relevant_var):

            # Saving directories
            counter = [self.counter_spec, self.counter_reac]

            save_path_max = self.reddir + '/MaxErrorReduction' + reduction_type + 'S' + str(counter[0]) + 'R' + \
                            str(counter[1])
            save_path_av = self.reddir + '/AvErrorReduction' + reduction_type + 'S' + str(counter[0]) + 'R' + \
                           str(counter[1])

            # Legend
            legend = []
            for index_x, x in enumerate(error_values):

                reactor_type_name = kwdict.reactor_labels[x[1]][0]

                if callable(x[2]):
                    error_name = x[2].__name__
                elif x[2] in reduced_mechanism.ctmech.species_names:
                    error_name = x[2]
                else:
                    error_name = kwdict.get_full_name(x[2])

                method_name = x[3]
                tolerence_string = x[4]

                separator = ', '

                if method_name == 'none':
                    method_name = ''
                    separator = ''

                if error_name == 'Velocity' and method_name == 'init':
                    error_name = 'Laminar flame speed'
                    method_name = ''
                    separator = ''

                legend.append(reactor_type_name + ' ' + error_name + separator + method_name
                              + ': ' + tolerence_string + '%')

            # Graph for maximum error
            backend = matplotlib.get_backend()
            matplotlib.use('Agg')
            title = 'Maximum error for ' + reduction_type + ' reduction'
            c = iter(cm.jet(np.linspace(0, 1, len(max_tol))))
            for id_y, y in enumerate(ydata_max):
                colour = next(c)
                graph.quick_graph(self.xdata, [yi * 100 for yi in y], xlabel=self.relevant_var, ylabel='error [%]',
                                  marker='x', line_style='-', colour=colour, title=title,
                                  legend=legend[id_y], reverse=True,
                                  show=False, save=None, hold=True)

                if self.max_error_plot:
                    graph.quick_graph(self.xdata,
                                      max_tol[id_y] * 100 * np.ones(len(self.xdata)),
                                      xlabel=self.relevant_var, ylabel='error [%]',
                                      ylim=2 * np.max(max_error_tol) * 100,
                                      marker=' ', line_style='--', colour=colour,
                                      title=title, reverse=True, show=False,
                                      save=None, hold=True)
            graph.arrow(getattr(reduced_mechanism, self.relevant_var),
                        np.max(max_tol) * 1.75 * 100, 0, -np.max(max_error_tol) * 0.5 * 100,
                        head_width=(np.max(self.xdata) - np.min(self.xdata)) / 100,
                        head_length=np.max(max_error_tol) * 0.1 * 100,
                        save=save_path_max)

            # Graph for average error
            title = 'Average error for ' + reduction_type + ' reduction'
            c = iter(cm.jet(np.linspace(0, 1, len(max_tol))))
            for id_y, y in enumerate(ydata_av):
                colour = next(c)
                graph.quick_graph(self.xdata, [yi * 100 for yi in y], xlabel=self.relevant_var, ylabel='error [%]',
                                  marker='x', line_style='-', colour=colour, title=title,
                                  legend=legend[id_y], reverse=True,
                                  show=False, save=None, hold=True)
                if self.max_error_plot:
                    graph.quick_graph(self.xdata,
                                      max_tol[id_y] * 100 * np.ones(len(self.xdata)),
                                      xlabel=self.relevant_var, ylabel='error [%]',
                                      ylim=2 * np.max(max_error_tol) * 100,
                                      marker=' ', line_style='--',  colour=colour,
                                      title=title, reverse=True, show=False,
                                      save=None, hold=True)

            graph.arrow(getattr(reduced_mechanism, self.relevant_var), np.max(max_error_tol) * 1.75 * 100, 0,
                        -np.max(max_error_tol) * 0.5 * 100,
                        head_width=(np.max(self.xdata) - np.min(self.xdata)) / 100,
                        head_length=np.max(max_error_tol) * 0.1 * 100,
                        save=save_path_av)

            logger.info('The maximal error against the number of ' + reduction_type
                        + ' has been saved as ' + save_path_max)
            logger.info('The average error against the number of ' + reduction_type
                        + ' has been saved as ' + save_path_av)
            logger.info('\n')
            matplotlib.use(backend)

    def simple_reduction(self, skeletal=False, sensitivity=0, permissive=False):
        """
        Linear reduction following species-reactions-QSS pattern

        Parameters
        ----------
        skeletal : bool
            if True does not do the QSS step (Default value = False)
        sensitivity : bool
            if True, performs the drgep-asa (aided sensitivity analysis)
        permissive : bool
            if True, continues computations even if one case of the list fails (Default value = False)

        Returns
        -------
        mechanism_list: list
            list of :func:`~ARCANE.mechanisms.Mechanism` for each reduction step (Default value = 0)
        """

        cases_list = self.cases_list
        mechanism = self.mechanism

        # Run simulations
        self.cantera_run_time += cases.run_cases(cases_list, mechanism, permissive=permissive, timing=True)

        # Checks for a super root mechanism as reference
        if not self.super_root_mechanism:
            self.super_root_mechanism = self.mechanism

        # Species reduction
        self.set_mechanism(mechanism)
        order = 'S' + ' : ' + str(mechanism.ns) + ' -> '
        species_reduced_mechanism = self.reduction('S', sensitivity=sensitivity, permissive=permissive)
        order = order + str(species_reduced_mechanism.ns) + '\n'

        # Reactions reduction
        self.set_mechanism(species_reduced_mechanism)
        order = order + 'R' + ' : ' + str(species_reduced_mechanism.nr) + ' -> '
        reactions_reduced_mechanism = self.reduction('R', sensitivity=sensitivity, permissive=permissive)
        order = order + str(reactions_reduced_mechanism.nr) + '\n'

        # QSS reduction
        if not skeletal:
            self.set_mechanism(reactions_reduced_mechanism)
            order = order + 'Q' + ' : ' + str(reactions_reduced_mechanism.nqss) + ' -> '
            qss_reduced_mechanism = self.reduction('QSS', permissive=permissive)
            order = order + str(qss_reduced_mechanism.nqss)
            mechanism_list = [mechanism, species_reduced_mechanism, reactions_reduced_mechanism, qss_reduced_mechanism]
        else:
            mechanism_list = [mechanism, species_reduced_mechanism, reactions_reduced_mechanism]

        logger.info('Here is the way your mechanism has been reduced :\n' + order)

        return mechanism_list

    def full_reduction(self, method='classic', full_mech_list=False, skeletal=False,
                       lumping=True, sensitivity=0, permissive=False):
        """
        Automatic full reduction with automatic method switch

        Parameters
        ----------
        method : str
            if 'classic', follows classical reduction, if 'layer', compute drgep by layer (Default value = 'classic')
        full_mech_list : bool
            if True outputs all the intermediate mechanisms (Default value = False)
        skeletal : bool
            if True, does not do the QSS step (Default value = False)
        lumping : bool
            if True, performs the lumping step (Default value = True)
        sensitivity : bool
            if True, performs the drgep-asa (aided sensitivity
            analysis) (Default value = 0)
        permissive : bool
            if True, continues computations even if one case of the list fails

        Returns
        -------
        mech_list : list
            list of class :func:`~ARCANE.mechanisms.Mechanism` objects for each reduction step (Default value = False)

        """

        # List of cases
        cases_list = self.cases_list

        # root_mechanism
        root_mechanism = self.mechanism

        # Setting empirically optimized error factors
        if self.error_factors == [0, 0, 0, 0]:
            self.error_factors = [0.8, 0.6, 0.8, 1]

            if skeletal:
                self.error_factors = [1, 0.8, 1, 1]

        # Checks for a super root mechanism as reference
        if not self.super_root_mechanism:
            self.super_root_mechanism = self.mechanism

        super_root_mechanism = self.super_root_mechanism

        # Run simulations for super root mechanism
        self.cantera_run_time += cases.run_cases(cases_list, super_root_mechanism, restore=self.restore,
                                                 permissive=permissive, timing=True)

        if root_mechanism != super_root_mechanism:
            # Run simulations for root mechanism
            self.cantera_run_time += cases.run_cases(cases_list, root_mechanism, restore=self.restore,
                                                     permissive=permissive, timing=True)

        # List of output mechanisms
        mech_list = [root_mechanism]
        mech_list_full = [root_mechanism]

        # Reference interest values
        root_ns = root_mechanism.ns
        root_nr = root_mechanism.nr
        root_nqss = root_mechanism.nqss

        # Storing the path
        reduction_type = ['None']
        reduction_path_ns = [root_ns]
        reduction_path_nr = [root_nr]
        reduction_path_nqss = [root_nqss]

        # Initializes the string telling the reduction order
        order = ''

        # Initial reduction being a species reduction
        if method == 'classic':
            disc = 2
            red_type = 'S'
            self.set_mechanism(root_mechanism)
            new_mechanism = self.reduction(red_type, sensitivity=sensitivity, permissive=permissive)
            self.counter_spec = self.counter_spec + 1
            order = 'S1 : ' + str(root_mechanism.ns) + ' -> ' + str(new_mechanism.ns) + '\n'

            mech_list_full.append(new_mechanism)

            # New relevant variables
            new_ns = new_mechanism.ns
            new_nr = new_mechanism.nr
            new_nqss = new_mechanism.nqss

        elif method == 'layer':
            disc = 6
            red_type = 'S'
            self.set_mechanism(root_mechanism)
            new_mechanism = root_mechanism
            # New relevant variables
            new_ns = new_mechanism.ns + 1
            new_nr = new_mechanism.nr + 1
            new_nqss = new_mechanism.nqss + 1

        else:
            new_ns = 0
            new_nr = 0
            new_nqss = 0
            disc = 0
            red_type = ''
            new_mechanism = ''
            logger.error('The method you are trying to test is not implemented !')

        if new_ns != reduction_path_ns[-1] \
                and new_nr != reduction_path_nr[-1] \
                and new_nqss != reduction_path_nqss[-1]:
            reduction_type.append(red_type)
            reduction_path_ns.append(new_ns)
            reduction_path_nr.append(new_nr)
            reduction_path_nqss.append(new_nqss)

        # If True, the reduction type as been switched from species to reduction or vice-versa
        switched = False

        error_dic_ref = []
        for case in cases_list:
            error_dic_ref.append(case.error_dict.copy())

        # Loops until all methods in every order have been tried
        for n in range(1, disc, 1):
            if disc > 2:
                logger.info("Here is error multiplied by :" + str(n / (disc - 1)))

            for id_case, case in enumerate(cases_list):
                for err_method in case.error_dict:
                    # Store the error method for the species variable
                    if err_method in kwdict.names['Y']:
                        for target in root_mechanism.targets:
                            if target != 'HeatRelease':
                                case.error_dict[err_method][1] = error_dic_ref[id_case][err_method][1] * n / (disc - 1)

                    # Store the error for the global variable
                    else:
                        case.error_dict[err_method] = error_dic_ref[id_case][err_method] * n / (disc - 1)

            # Force to do the first loop
            if n > 1:
                new_ns = 0
                new_nr = 0

            while True:

                # If the number of species changed, reduce species again
                if root_ns != new_ns:

                    red_type = 'S'

                    if self.step == 'auto':
                        if len(self.ns_auto_stepping) > 1:
                            self.ns_auto_stepping = self.ns_auto_stepping[1:]
                        else:
                            self.step = 1

                    switched = False

                # If the number of reactions changed and not the number of species, reduce reactions again
                elif root_ns == new_ns and root_nr != new_nr:

                    red_type = 'R'

                    if self.step == 'auto':
                        if len(self.nr_auto_stepping) > 1:
                            self.nr_auto_stepping = self.nr_auto_stepping[1:]
                        else:
                            self.step = 1

                    switched = False

                # If the same reduction gave consecutively the same result, switch method
                else:
                    if self.step == 'auto':
                        if red_type == 'S':
                            if len(self.ns_auto_stepping) > 1:
                                self.ns_auto_stepping = self.ns_auto_stepping[1:]
                                switched = False
                            else:
                                self.step = 1
                        else:
                            if len(self.nr_auto_stepping) > 1:
                                self.nr_auto_stepping = self.nr_auto_stepping[1:]
                                switched = False
                            else:
                                self.step = 1

                    # If method switched at previous iteration, ends the reduction
                    if switched:
                        break
                    # Switching the method
                    else:
                        #  Switching the reduction ype
                        if red_type == 'S':
                            red_type = 'R'
                        else:
                            red_type = 'S'

                        switched = True

                if red_type == 'S':
                    self.counter_spec = self.counter_spec + 1
                    order = order + 'S' + str(self.counter_spec) + ' : '
                else:
                    self.counter_reac = self.counter_reac + 1
                    order = order + 'R' + str(self.counter_reac) + ' : '

                # Previous mechanism becomes the reference for change assessing
                root_ns = new_ns
                root_nr = new_nr

                # Setting and launching new reduction
                self.set_mechanism(new_mechanism)
                if red_type == 'S':
                    order = order + str(new_mechanism.ns) + ' -> '
                else:
                    order = order + str(new_mechanism.nr) + ' -> '

                new_mechanism = self.reduction(red_type, sensitivity=sensitivity, permissive=permissive)

                if red_type == 'S':
                    order = order + str(new_mechanism.ns) + '\n'
                else:
                    order = order + str(new_mechanism.nr) + '\n'

                new_ns = new_mechanism.ns
                new_nr = new_mechanism.nr
                new_nqss = new_mechanism.nqss

                if new_ns != reduction_path_ns[-1] \
                    or new_nr != reduction_path_nr[-1] \
                    or new_nqss != reduction_path_nqss[-1]:
                    reduction_type.append(red_type)
                    reduction_path_ns.append(new_ns)
                    reduction_path_nr.append(new_nr)
                    reduction_path_nqss.append(new_nqss)

                mech_list_full.append(new_mechanism)

            mech_list.append(new_mechanism)

        if lumping:
            # Previous mechanism becomes the reference for change assessing
            root_ns = new_ns

            #  Lumping reduction
            red_type = 'L'
            order = order + red_type + ' : ' + str(new_mechanism.ns) + ' -> '
            self.set_mechanism(new_mechanism)

            new_mechanism = self.lumping_reduction(permissive=permissive)

            order = order + str(new_mechanism.ns) + '\n'
            new_ns = new_mechanism.ns
            new_nr = new_mechanism.nr
            new_nqss = new_mechanism.nqss

            if new_ns != reduction_path_ns[-1] \
                    or new_nr != reduction_path_nr[-1] \
                    or new_nqss != reduction_path_nqss[-1]:
                reduction_type.append(red_type)
                reduction_path_ns.append(new_ns)
                reduction_path_nr.append(new_nr)
                reduction_path_nqss.append(new_nqss)

            mech_list_full.append(new_mechanism)
            mech_list.append(new_mechanism)

            # DRGEP loop
            if root_ns != new_ns:
                # Initial reduction being a species reduction
                red_type = 'S'
                self.set_mechanism(new_mechanism)
                self.counter_spec = self.counter_spec + 1
                order = order + 'S' + str(self.counter_spec) + ' : '
                order = order + str(new_mechanism.ns) + ' -> '

                new_mechanism = self.reduction(red_type, sensitivity=sensitivity, permissive=permissive)

                order = order + str(new_mechanism.ns) + '\n'

                # New relevant variables
                new_ns = new_mechanism.ns
                new_nr = new_mechanism.nr

                mech_list_full.append(new_mechanism)

                # If True, the reduction type as been switched from species to reduction or vice-versa
                switched = False

                # Loops until all methods in every order has been tried
                while True:
                    # If the number of species changed, reduce species again
                    if root_ns != new_ns:

                        red_type = 'S'
                        switched = False

                    # If the number of reactions changed and not the number of species, reduce reactions again
                    elif root_ns == new_ns and root_nr != new_nr:

                        red_type = 'R'
                        switched = False

                    # If the same reduction gave consecutively the same result, switch method
                    else:
                        # If method switched at previous iteration, ends the reduction
                        if switched:
                            break
                        # Switching the method
                        else:
                            #  Switching the reduction ype
                            if red_type == 'S':
                                red_type = 'R'
                            else:
                                red_type = 'S'

                            switched = True

                    if red_type == 'S':
                        self.counter_spec = self.counter_spec + 1
                        order = order + 'S' + str(self.counter_spec) + ' : '
                    else:
                        self.counter_reac = self.counter_reac + 1
                        order = order + 'R' + str(self.counter_reac) + ' : '

                    # Previous mechanism becomes the reference for change assessing
                    root_ns = new_ns
                    root_nr = new_nr

                    # Setting and launching new reduction
                    self.set_mechanism(new_mechanism)
                    if red_type == 'S':
                        order = order + str(new_mechanism.ns) + ' -> '
                    else:
                        order = order + str(new_mechanism.nr) + ' -> '

                    new_mechanism = self.reduction(red_type, sensitivity=sensitivity, permissive=permissive)

                    if red_type == 'S':
                        order = order + str(new_mechanism.ns) + '\n'
                    else:
                        order = order + str(new_mechanism.nr) + '\n'

                    new_ns = new_mechanism.ns
                    new_nr = new_mechanism.nr
                    new_nqss = new_mechanism.nqss

                    if new_ns != reduction_path_ns[-1] \
                    or new_nr != reduction_path_nr[-1] \
                    or new_nqss != reduction_path_nqss[-1]:
                        reduction_type.append(red_type)
                        reduction_path_ns.append(new_ns)
                        reduction_path_nr.append(new_nr)
                        reduction_path_nqss.append(new_nqss)

                    mech_list_full.append(new_mechanism)

                mech_list.append(new_mechanism)

        if not skeletal:
            #  QSS reduction
            red_type = 'QSS'
            order = order + 'QSS : ' + str(new_mechanism.nqss) + ' -> '
            self.set_mechanism(new_mechanism)
            new_mechanism = self.reduction(red_type, permissive=permissive)
            order = order + str(new_mechanism.nqss)

            new_ns = new_mechanism.ns
            new_nr = new_mechanism.nr
            new_nqss = new_mechanism.nqss

            if new_ns != reduction_path_ns[-1] \
                    or new_nr != reduction_path_nr[-1] \
                    or new_nqss != reduction_path_nqss[-1]:
                reduction_type.append(red_type)
                reduction_path_ns.append(new_ns)
                reduction_path_nr.append(new_nr)
                reduction_path_nqss.append(new_nqss)

            mech_list_full.append(new_mechanism)
            mech_list.append(new_mechanism)

        if full_mech_list:
            mech_list = mech_list_full

        # Timing the reduction
        final_time = time.time()
        elapsed_time = final_time - self.start_time
        elapsed_time_string = tools.format_time(elapsed_time)

        cantera_total_run_time = sum(self.cantera_run_time)

        cantera_run_time_string = tools.format_time(cantera_total_run_time)
        red_run_time_string = tools.format_time(elapsed_time - cantera_total_run_time)

        run_time_percentage = cantera_total_run_time * 100 / elapsed_time
        run_time_percentage_string = str(int(run_time_percentage))
        red_time_percentage_string = str(100 - int(run_time_percentage))

        logger.info('The reduction took ' + elapsed_time_string)
        logger.info('Time spent running cases: ' + cantera_run_time_string
                    + ' (' + run_time_percentage_string + '%)')
        logger.info('Time spent in the reduction algorithm: '
                    + red_run_time_string + ' (' + red_time_percentage_string + '%)\n')

        logger.info('The mechanism has been reduced in the following order:\n' + order)

        reduction_dict = {'type': reduction_type,
                          'ns': reduction_path_ns,
                          'nr': reduction_path_nr,
                          'nqss': reduction_path_nqss}
        logger.info('\n\nHere is the reduction path dictionary:\n' + str(reduction_dict))

        return mech_list

def get_head_and_tail(reduction_type):
    """Get the head and tail for the reduction type chosen

    Parameters
    ----------
    reduction_type : str
        reduction type can be either 'species', 'reactions', 'qss', 'lumping'

    Returns
    -------
    head : str
        head of display
    tail : str
        tail of display
    """
    if reduction_type in ['species', 'S']:
        head, tail = display.head_tail_charac(text='Species reduction', size='big', style='#')
    elif reduction_type in ['reactions', 'R']:
        head, tail = display.head_tail_charac(text='Reactions reduction', size='big', style='#')
    elif reduction_type in ['QSS', 'qss']:
        head, tail = display.head_tail_charac(text='QSS reduction', size='big', style='#')
    elif reduction_type in ['lumping', 'L']:
        # CHECK IF OK
        head, tail = display.head_tail_charac(text='Species lumping', size='big', style='#')
    else:
        head = 0
        tail = 0
        logger.error('ERROR: Reduction type is not valid ! Valid ones are : species, reactions and QSS')
    logger.info(head)
    return head, tail
