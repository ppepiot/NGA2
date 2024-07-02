# Module for generating a pyrolysis mechanism
import matplotlib.pyplot as plt

import ARCANE.cases as cases
import ARCANE.display as display
import ARCANE.sampling as sampling
import ARCANE.postproc as postproc
import ARCANE.tools as tools
import ARCANE.display as display

import sys
import scipy.optimize as optimize
import numpy as np
import warnings
import scipy.optimize as opt

np.set_printoptions(threshold=sys.maxsize)
warnings.filterwarnings('ignore')

logger = display.Logger()
logger.set_log('logPathway')
logger.create_log('logPathway')


class Pyrolysis:

    def __init__(self, cases_list, mechanism, exc_species_list):
        """Initialisation of the Pathway object that will contain everything needed for the lumping

        Parameters
        ----------
        cases_list : list
            list of class :func:`~ARCANE.cases.Case` objects
        mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            class :func:`~ARCANE.mechanisms.Mechanism` object to reduce

        """

        # Initialisation of important parameters
        self.mechanism = mechanism

        # Reactions that contain each of the species from exc_species_list should not be lumped
        self.exc_species_list = exc_species_list

        # Separating the cases for reduction from the complete list
        self.all_cases_list = cases_list
        self.cases_list = [case for case in cases_list if '0D' in case.reactor_type]

        # Running the cases
        cases.run_cases(cases_list, mechanism)

        # Initial number of reactions
        self.nr = self.mechanism.nr

        # Number of samples for each case
        self.n_samples = 200

        # Sets if the Arrhenius fit is done with the temperature correction
        self.temperature_correction_Arrhenius = False

        # List of final lumped reactions
        self.lumped_reactions_objects = []

        # Retrieving fuel species
        self.fuel_species_dict = self.cases_list[0].fuel
        logger.info(f"\nPathway lumping performed for fuel: {self.fuel_species_dict}")
        one_composition = True
        for case in self.cases_list:
            if case.fuel != self.fuel_species_dict:
                one_composition = False
                logger.warning(f"Different fuel composition {case.fuel} found in the cases list, "
                               f"this will lead to errors.")

        if not one_composition:
            quit()

        self.fuel_species = list(self.fuel_species_dict.keys())

        # Default importance thresholds
        self.reactants_importance_threshold = [1e-3] * len(self.fuel_species)
        self.products_importance_threshold = [1e-3] * len(self.fuel_species)

        # Default number of reactions and products
        self.number_of_reactions = [1] * len(self.fuel_species)
        self.number_of_products = [4] * len(self.fuel_species)

        # Choosing if the reactions and products selection will be based on thresholds or number
        self.use_threshold_for_reactions = True
        self.use_threshold_for_products = True

        self.species_objects = self.mechanism.ctmech.species
        # Creating a composition dictionary
        composition_dict = {}
        for spec in self.mechanism.species_names:
            composition_dict[spec] = {'C': self.mechanism.ctmech.n_atoms(spec, 'C'),
                                      'H': self.mechanism.ctmech.n_atoms(spec, 'H'),
                                      'O': self.mechanism.ctmech.n_atoms(spec, 'O')}

        # Surrogate composition
        composition_dict['Surrogate'] = {'C': 0, 'H': 0, 'O': 0}
        for fuel in self.fuel_species_dict:
            composition_dict['Surrogate']['C'] += composition_dict[fuel]['C'] * self.fuel_species_dict[fuel]
            composition_dict['Surrogate']['H'] += composition_dict[fuel]['H'] * self.fuel_species_dict[fuel]
            composition_dict['Surrogate']['O'] += composition_dict[fuel]['O'] * self.fuel_species_dict[fuel]

        self.composition_dict = composition_dict

        self.discard_species_in_build = True

        self.must_be_irreversible = True

        # If true, one of the reactions will be the fuel homolysis
        self.fuel_homolysis = True

    def set_reactants_importance_threshold(self, value):
        """Setting the list of importance thresholds for the H abstraction species

        Parameters
        ----------
        value :
            value of the threshold either in float or in list of multi-fuel is used

        """
        if type(value) != list:
            value = [value] * len(self.fuel_species)

        self.reactants_importance_threshold = value

    def set_products_importance_threshold(self, value):
        """Setting the list of importance thresholds for the products

        Parameters
        ----------
        value :
            value of the threshold either in float or in list of multi-fuel is used

        """
        if type(value) != list:
            value = [value] * len(self.fuel_species)

        self.products_importance_threshold = value

    def get_pyrolysis_species(self):
        """The species and reactions are sorted according to their role in the mechanism

        Returns
        -------
        self.reactants : list
            list of species identified as reactants for the pyrolysis lumping
        self.products : list
            list of species identified as products for the pyrolysis lumping
        self.precursors : list
            list of species identified as reactants but below the threshold
        self.byproducts : list
            list of species identified as products but below the threshold
        self.reactions : list
            list of reactions from self.mechanism.ctmech.reactions() (class `Cantera.Solution.reactions`) identified as lumpable
        self.reactions_indices : list
            list of indices for the reactions to be lumped


        """

        # Identifying the species that are part of the pyrolysis process
        # These species are defined as fragments of the fuel species
        # They come from irreversible reactions where the fuel or a fuel fragment is broken down
        # The species will be selected as the biggest species being a product of a pyrolysis reaction
        # With a carbon atom greater than 4

        # In order to avoid several loops, all the species with more than 4 carbon atoms are selected a priori
        # The ones not following the previous rules will be discarded
        species_objects = self.mechanism.ctmech.species
        potential_pyrolysis_species = [spec.name for spec in species_objects() if 'C' in spec.composition
                                       and spec.composition['C'] > 4]

        final_pyrolysis_reactants = []
        final_pyrolysis_products = []
        pyrolysis_precursors = []
        pyrolysis_byproducts = []

        # Saving the pyrolysis reaction (with and without the extra species)
        pyrolysis_reactions = []
        pyrolysis_reactions_expressions = []
        pyrolysis_reactions_indices = []

        reaction_objects = self.mechanism.ctmech.reactions()

        def dummy_function(float):
            """

            Parameters
            ----------
            float :
                dummy argument


            Returns
            -------
            False : bool
                dummy return False

            """
            return False

        if self.must_be_irreversible:
            is_reversible_function = self.mechanism.ctmech.is_reversible
        else:
            is_reversible_function = dummy_function

        condition_for_reactant = "(species_objects(reactant).composition['C'] > species_objects(" \
                    "pyrolysis_reactant).composition['C'] or (species_objects(reactant).composition['C'] == " \
                    "species_objects(pyrolysis_reactant).composition['C'] and species_objects(" \
                    "reactant).composition['H'] > species_objects(pyrolysis_reactant).composition['H'])) "

        #condition_for_reactant = "species_objects(reactant).composition['C'] > species_objects(" \
        #            "pyrolysis_reactant).composition['C'] "

        # Verification of the species
        for index_reac, reac in enumerate(reaction_objects):
            reactants = reac.reactants
            reactants_list = list(reactants.keys())
            products = reac.products
            products_list = list(products)

            # Reaction must not contain any of the exceptional species
            if len(set(reactants_list + products_list).intersection(set(self.exc_species_list))) == 0:
                # Reaction must be irreversible and have a potential species as reactant
                reactant_in_potential_species = list(set(reactants_list) & set(potential_pyrolysis_species))
                if not is_reversible_function(index_reac) and reactant_in_potential_species:
                    # Retrieving the biggest species among the reactants
                    pyrolysis_reactant = reactant_in_potential_species[0]
                    # If several pyrolysis species are present, take the biggest
                    if len(reactant_in_potential_species) > 1:
                        for reactant in reactant_in_potential_species:
                            if eval(condition_for_reactant):
                                pyrolysis_reactant = reactant

                    # Determining if there is a pyrolysis product
                    pyrolysis_products = []

                    condition_for_product = "(species_objects(prod).composition['C'] < species_objects(" \
                                "pyrolysis_reactant).composition['C'] or (species_objects(prod).composition['C'] == " \
                                "species_objects(pyrolysis_reactant).composition['C'] and species_objects(" \
                                "prod).composition['H'] < species_objects(pyrolysis_reactant).composition['H'])) "

                    #condition_for_product = "species_objects(prod).composition['C'] < species_objects(" \
                    #            "pyrolysis_reactant).composition['C'] "

                    for prod in products_list:
                        if 'C' in species_objects(prod).composition and eval(condition_for_product):
                            pyrolysis_products.append(prod)

                    if pyrolysis_products:
                        # Storing precursors and byproducts of pyrolysis
                        pyrolysis_precursors += reactants_list
                        pyrolysis_byproducts += products_list

                        pyrolysis_product = pyrolysis_products[0]
                        # If there is several fitting products, takes the biggest
                        if len(pyrolysis_products) > 2:
                            pyrolysis_product = pyrolysis_products[0]
                            for spec in pyrolysis_products[1:]:
                                if species_objects(spec).composition['C'] > \
                                        species_objects(pyrolysis_product).composition['C']:
                                    pyrolysis_product = spec

                        # Saving the results
                        final_pyrolysis_reactants.append(pyrolysis_reactant)
                        final_pyrolysis_products.append(pyrolysis_product)
                        # Reactions objects
                        pyrolysis_reactions.append(reac)
                        # Reactions indices
                        pyrolysis_reactions_indices.append(index_reac)
                        # Expressions
                        pyrolysis_reactions_expressions.append(reac.equation)

        # Removing doubles
        final_pyrolysis_reactants = list(set(final_pyrolysis_reactants))
        final_pyrolysis_products = list(set(final_pyrolysis_products) - set(final_pyrolysis_reactants))
        pyrolysis_precursors = list(set(pyrolysis_precursors) - set(final_pyrolysis_reactants))
        pyrolysis_byproducts = list(
                set(pyrolysis_byproducts) - set(final_pyrolysis_products) - set(final_pyrolysis_reactants))

        # Storing the results
        self.reactants = final_pyrolysis_reactants
        self.products = final_pyrolysis_products
        self.precursors = pyrolysis_precursors
        self.byproducts = pyrolysis_byproducts
        self.reactions = pyrolysis_reactions
        self.reactions_expressions = pyrolysis_reactions_expressions
        self.reactions_indices = pyrolysis_reactions_indices

        # print(f"Pyrolysis reactants: {final_pyrolysis_reactants}")
        # print(f"Pyrolysis products: {final_pyrolysis_products}")
        # print(f"Pyrolysis precursors: {pyrolysis_precursors}")
        # print(f"Pyrolysis byproducts: {pyrolysis_byproducts}")
        # print(f"Pyrolysis reactions:")
        # for reac in pyrolysis_reactions_expressions:
        #     print(reac)
        # quit()

        return self.reactants, self.products, self.precursors, self.byproducts, self.reactions, self.reactions_indices

    def load_samples(self):
        """Loading the samples database

        """

        self.samples_database = sampling.create_samples_database(self.cases_list, self.mechanism,
                                                                 n_samples=self.n_samples, time_normalised=True,
                                                                 integrated=False)

        self.samples_database = postproc.extend_data_dict(self.samples_database)

        self.n_samples = len(self.samples_database['grid'])

        self.fuel_concentration = np.zeros(self.n_samples)
        for fuel in self.fuel_species_dict:
            self.fuel_concentration += self.samples_database[f"c_{fuel}"]

        self.samples_database[f'c_Surrogate'] = self.fuel_concentration

    def sorted_source_terms(self):
        """Retrieves the data from the computations and sort them into consumption and creation of species
        in pyrolysis and oxidation reactions

        """

        reactions = self.mechanism.reactions
        species = self.mechanism.species_names

        pyrolysis_reactions = self.reactions_indices

        pyro_con = {}
        pyro_prod = {}
        pyro_net = {}
        oxi_con = {}
        oxi_prod = {}

        # Extending the data_dict
        reaction_rates = self.samples_database['net rates of progress']

        for spec in species:
            pyro_con[spec] = np.zeros([self.n_samples])
            pyro_prod[spec] = np.zeros([self.n_samples])
            oxi_con[spec] = np.zeros([self.n_samples])
            oxi_prod[spec] = np.zeros([self.n_samples])

        # Rates separation
        for i_reac, reac in enumerate(reactions):
            reactants = list(reac.reactants.keys())
            products = list(reac.products.keys())
            for i_spec, spec in enumerate(species):
                if spec in reactants + products:
                    # Determining if pyrolysis or not
                    if i_reac in pyrolysis_reactions:
                        # If in reactants: consumed
                        if spec in reactants:
                            pyro_con[spec] += reaction_rates[i_reac, :]
                        # If in products: produced
                        else:
                            pyro_prod[spec] += reaction_rates[i_reac, :]
                    # Oxidation reactions
                    else:
                        # If in reactants: consumed
                        if spec in reactants:
                            oxi_con[spec] += reaction_rates[i_reac, :]
                        # If in products: produced
                        else:
                            oxi_prod[spec] += reaction_rates[i_reac, :]

        # Cracking consumption rates
        self.cracking_net_rates_dict = pyro_net
        self.cracking_consumption_rates_dict = pyro_con
        self.cracking_production_rates_dict = pyro_prod

    def compute_hydrogen_abstraction_importance(self):
        """Computes the alpha coefficient for each precursor


        TODO
        ----
        Put a warning for when the alpha values exceeds 1 telling that the case is not representative


        """

        alpha_dict = {}
        self.normalised_time = self.samples_database['normalised grid interval']

        fuel_decomposition_sum = 0
        for fuel in self.fuel_species_dict:
            fuel_decomposition = np.array(self.cracking_consumption_rates_dict[fuel]) * self.normalised_time
            fuel_decomposition_sum += sum(fuel_decomposition)

        for spec in self.precursors:
            spec_abstraction = np.array(self.cracking_consumption_rates_dict[spec]) * self.normalised_time
            spec_abstraction_sum = sum(spec_abstraction)
            if fuel_decomposition_sum > 0:
                alpha = spec_abstraction_sum / fuel_decomposition_sum
            else:
                alpha = 0

            if alpha > 1e-05:
                alpha_dict[spec] = alpha

        self.precursors_coefficients = alpha_dict

    def compute_hydrogen_abstraction_rates(self):
        """Computes consumption rates of the identified precursors

        """
        import matplotlib.pyplot as plt

        hydro_abst_rate_coeffs = {}
        sum_hydro_abst_rates = np.zeros(self.n_samples)

        for index_spec, spec in enumerate(self.selected_precursors):
            hydro_abst_rate_coeffs[spec] = self.cracking_consumption_rates_dict[spec] \
                                           / (self.samples_database[f"c_{spec}"] * self.fuel_concentration)

            # Discarding values arising from numerical errors
            for fuel in self.fuel_species_dict:
                hydro_abst_rate_coeffs[spec][self.samples_database[f"c_{fuel}"] < 1e-10] = np.nan

            sum_hydro_abst_rates += self.cracking_consumption_rates_dict[spec]

        self.sum_hydro_abst_rates = sum_hydro_abst_rates

        self.hydrogen_abstraction_rate_coefficients = hydro_abst_rate_coeffs
        # plt.plot(1000/np.array(self.samples_database[f"T"]), np.log(hydro_abst_rates[spec]), 'o', label=spec)

        # plt.yscale('log')
        # plt.legend()
        # plt.show()

    def compute_fuel_decomposition_rate(self):
        """Computes the consumption rate of the fuel as
        k(T) = cdot(T) / c_fuel(T)

        """

        fuel_consumption_rate = np.zeros(self.n_samples)
        for fuel in self.fuel_species_dict:
            index_fuel = self.mechanism.species_names.index(fuel)
            fuel_consumption_rate += - self.samples_database[f"net production rates"][index_fuel, :]
            # fuel_consumption_rate += self.cracking_consumption_rates_dict[fuel]

        modified_fuel_decomposition_rate = fuel_consumption_rate - self.sum_hydro_abst_rates

        rate_coefficient = modified_fuel_decomposition_rate / self.fuel_concentration

        # Discarding zones where the fuel concentration is too low
        for fuel in self.fuel_species_dict:
            rate_coefficient[self.samples_database[f"c_{fuel}"] < 1e-5] = np.nan
        # Discarding negative values (arising from faulty reactions identification)
        rate_coefficient[rate_coefficient < 0] = np.nan

        self.fuel_decomposition_rate_coefficient = rate_coefficient

        # plt.plot(self.samples_database[f"T"], modified_fuel_decomposition_rate, 'o')
        # index_fuel = self.mechanism.species_names.index(self.fuel_species)
        # plt.plot(self.samples_database[f"T"], - self.samples_database[f"net production rates"][index_fuel, :], 'o')
        # plt.show()
        # plt.plot(self.samples_database[f"T"], self.samples_database[f"c_{self.fuel_species}"], 'o')
        # plt.show()
        # plt.plot(self.samples_database[f"T"], rate_coefficient, 'o')
        # plt.show()
        # for spec in self.precursors:
        #     print(spec)
        #     new_rate = modified_fuel_decomposition_rate / self.samples_database[f"c_{spec}"]
        #     plt.plot(self.samples_database[f"T"], np.log(new_rate))
        #     plt.show()
        # for spec in self.precursors:
        #     print(spec, ' + XYLENE')
        #     new_rate = modified_fuel_decomposition_rate / ( self.samples_database[f"c_{self.fuel_species}"] * self.samples_database[f"c_{spec}"])
        #     plt.plot(self.samples_database[f"T"], np.log(new_rate))
        #     plt.show()
        # plt.plot(1000/np.array(self.samples_database[f"T"]), np.log(rate_coefficient), 'o')
        # plt.plot(1000 / np.array(self.samples_database[f"T"]), modified_fuel_decomposition_rate, 'o')
        # plt.show()
        # quit()

    def fit_rate_coeffcients(self):
        """Fit the reaction rates with an Arrhenius function

        """

        Arrhenius_parameters_dict = dict.fromkeys(self.rate_coefficients.keys(), [])

        for spec in Arrhenius_parameters_dict:
            rate_coefficients_untreated = self.rate_coefficients[spec]

            # Pre-treat arrays to discard NaN values
            temperatures = np.array(self.samples_database[f"T"])[~np.isnan(rate_coefficients_untreated)]
            rate_coeffs = rate_coefficients_untreated[~np.isnan(rate_coefficients_untreated)]
            # Pre-treat arrays to discard infs values
            temperatures = temperatures[~np.isinf(rate_coeffs)]
            rate_coeffs = rate_coeffs[~np.isinf(rate_coeffs)]
            # Pre-treat arrays to discard values not valid for log (< 0)
            temperatures = temperatures[rate_coeffs > 0]
            rate_coeffs = rate_coeffs[rate_coeffs > 0]

            # Fitting values
            # self.temperature_correction_Arrhenius = True
            if self.temperature_correction_Arrhenius:
                reverse_Arrh, pcov = opt.curve_fit(tools.calculate_lnk, temperatures, np.log(rate_coeffs))
            else:
                reverse_Arrh, pcov = opt.curve_fit(tools.calculate_lnk_simple, temperatures, np.log(rate_coeffs))

            # plt.plot(temperatures)
            # plt.plot(temperatures, np.log(rate_coeffs), 'o', label=spec)
            # if self.temperature_correction_Arrhenius:
            #     plt.plot(temperatures, tools.calculate_lnk(temperatures, *reverse_Arrh), label=spec + ' fit')
            # else:
            #     plt.plot(temperatures, tools.calculate_lnk_simple(temperatures, *reverse_Arrh), label=spec + ' fit')

            # # Checking for extreme values
            # if any(abs(reverse_Arrh) > 100):
            #     logger.error(f'The data cannot be fitted properly, {reverse_Arrh}')
            #     quit()

            # Retrieving real values
            reverse_Arrh[0] = np.exp(reverse_Arrh[0])
            if self.temperature_correction_Arrhenius:
                reverse_Arrh[2] = reverse_Arrh[2] * 8314.4621
            else:
                reverse_Arrh[1] = reverse_Arrh[1] * 8314.4621

            Arrhenius_parameters_dict[spec] = reverse_Arrh

        self.Arrhenius_parameters_dict = Arrhenius_parameters_dict

        # print(self.Arrhenius_parameters_dict)
        # plt.legend()
        # plt.show()
        # quit()

    def identify_products(self):
        """Identifies the potential products for the lumped reactions
        The identification is done by looking at the production of products and byproducts in pyrolysis reactions

        """

        # Find the peak rates of fuel consumption for each case in the sampled data
        list_of_cases = list(set(self.samples_database['case id']))
        local_max = []
        for case_id in list_of_cases:
            local_ids = [local_id == case_id for local_id in self.samples_database['case id']]
            for fuel in self.fuel_species_dict:
                local_max.append(list(self.cracking_consumption_rates_dict[fuel]).index(
                        max(self.cracking_consumption_rates_dict[fuel][local_ids])))
        local_max = list(set(local_max))

        # Computing the sum over all species of the production rates
        sum_of_productions = 0
        for spec in self.products + self.byproducts:
            for index in local_max:
                sum_of_productions += self.cracking_production_rates_dict[spec][index] * self.normalised_time[index]

        # Computing each species coefficient
        products_importance_coefficients = {}
        for spec in self.products:  # + self.byproducts:
            products_importance_coefficients[spec] = 0
            for index in local_max:
                products_importance_coefficients[spec] += self.cracking_production_rates_dict[spec][index] \
                                                          * self.normalised_time[index] / sum_of_productions

        self.products_coefficients = products_importance_coefficients

    def solve_products_coefficients(self):
        """Resolves the coupled linear systems of the products productions and the atomic conservation

        """

        # Parameters retrieval
        all_products = [value for sub_list in self.products_by_reaction.values() for value in sub_list]
        total_number_of_products = len(all_products)
        unique_products = list(set(all_products))
        number_of_unique_products = len(unique_products)

        # Creation of the first system
        # The sum of the reaction rates of each species across all reactions must be equal to its total production rate
        # This takes the form of a Ax = b problem
        # A holds the (weighted) reaction rates (at each sample for each unique species (nT*nP)
        # for each species in each reaction (nP(R)*nR))
        # b holds the (weighted) production rates (at each sample for each unique species (nT*nP))
        # Each quantity is weighted by the normalized time of the sample divided by
        # the mean concentration of the species (at the sample)

        # Only using relevant data points
        too_low_concentration_mask = self.fuel_concentration > 1e-10
        temperatures = np.array(self.samples_database[f"T"])[too_low_concentration_mask]
        local_n_samples = len(temperatures)

        # Initialising the matrices
        reaction_rates_matrix = np.zeros((local_n_samples * number_of_unique_products, total_number_of_products))
        production_rates_vector = np.zeros((local_n_samples * number_of_unique_products))
        stoichiometric_coefficients_vector = np.zeros(total_number_of_products)

        # Filling the matrices
        counter = 0
        for index_reaction, reactant in enumerate(self.products_by_reaction):

            # Computing reaction rate from fitted Arrhenius and detailed concentrations
            if self.temperature_correction_Arrhenius:
                rate_coefficient = self.Arrhenius_parameters_dict[reactant][0] \
                                   * temperatures ** self.Arrhenius_parameters_dict[reactant][1] \
                                   * np.exp(-self.Arrhenius_parameters_dict[reactant][2] / (8.3144621 * temperatures))
            else:
                rate_coefficient = self.Arrhenius_parameters_dict[reactant][0] \
                                   * np.exp(-self.Arrhenius_parameters_dict[reactant][1] / (8.3144621 * temperatures))
            reaction_rate = rate_coefficient

            for spec in self.reactants_by_reaction[reactant]:
                spec_concentration = self.samples_database[f'c_{spec}'][too_low_concentration_mask]
                reaction_rate *= spec_concentration

            # Loop on the products
            local_counter = 0
            for spec in self.products_by_reaction[reactant]:
                unique_species_index = unique_products.index(spec)

                start_index = unique_species_index * local_n_samples
                end_index = (unique_species_index + 1) * local_n_samples

                index_product_at_reaction = counter + local_counter

                reaction_rates_matrix[start_index:end_index, index_product_at_reaction] = reaction_rate

                local_counter += 1

            counter += local_counter

        for unique_species_index, spec in enumerate(unique_products):

            global_species_index = self.mechanism.species_names.index(spec)

            if spec not in self.selected_precursors:
                net_production_rates = self.samples_database['creation rates'][global_species_index,
                                                                               too_low_concentration_mask]
            else:
                net_production_rates = self.samples_database['net production rates'][global_species_index,
                                                                                     too_low_concentration_mask]

            start_index = unique_species_index * local_n_samples
            end_index = (unique_species_index + 1) * local_n_samples

            production_rates_vector[start_index:end_index] = net_production_rates

        # Creation of the second system
        # The sum of number of elements of each species for each reactions must be equal to
        # the total number of this element arsing for each species
        # This takes the form of a Cx = d problem
        # C holds the number of element in the reaction (in each reaction for each element in this reaction (nR*nE(R))
        # for each species in each reaction (nP(R)*nR))
        # d holds the total number of element arsing from each species (for each element
        # for each unique species (nE(R)*nP))

        elements_per_reaction = {}
        total_elements_per_reaction = {}
        for reaction in self.reactants_by_reaction:
            local_elements = []
            total_elements = {}
            for reactant in self.reactants_by_reaction[reaction]:
                for element in self.composition_dict[reactant]:
                    if self.composition_dict[reactant][element] > 0:
                        local_elements.append(element)
                        if element not in total_elements:
                            total_elements[element] = self.composition_dict[reactant][element]
                        else:
                            total_elements[element] += self.composition_dict[reactant][element]

            elements_per_reaction[reaction] = list(set(local_elements))
            total_elements_per_reaction[reaction] = total_elements

        number_of_elements_per_reaction = [len(elements_per_reaction[spec]) for spec in elements_per_reaction]

        total_number_of_elements = sum(number_of_elements_per_reaction)

        element_per_reaction_matrix = np.zeros((total_number_of_elements, total_number_of_products))
        element_per_species_vector = np.zeros((total_number_of_elements))

        counter_prod = 0
        counter_element = 0
        for index_reaction, reaction in enumerate(self.products_by_reaction):

            local_counter_element = 0
            for index_element, element in enumerate(elements_per_reaction[reaction]):

                index_element_in_reac = counter_element + local_counter_element

                local_counter_prod = 0
                for product in self.products_by_reaction[reaction]:
                    index_element_in_spec = counter_prod + local_counter_prod

                    element_per_reaction_matrix[index_element_in_reac,
                                                index_element_in_spec] = self.composition_dict[product][element]

                    local_counter_prod += 1

                element_per_species_vector[index_element_in_reac] = total_elements_per_reaction[reaction][element]

                local_counter_element += 1

            counter_element += local_counter_element
            counter_prod += local_counter_prod

        # Standard Least Squares
        def lsqmin(x, a, b):
            """Method to compute least squares

            Parameters
            ----------
            x :
                variable vector

            a :
                constant vector

            b :
                constant vector


            Returns
            -------
             : ndarray
                norme de l'expression A*X - B

            """
            return np.linalg.norm(np.matmul(a, x) - b)

        # --------- Call COBYLA routine to enforce constraint and obtain stoich. coeffs.
        cons = [{'type': 'ineq',
                 'fun' : lambda x: np.matmul(element_per_reaction_matrix, x) - np.array(element_per_species_vector)},
                {'type': 'ineq',
                 'fun' : lambda x: np.array(element_per_species_vector) - np.matmul(element_per_reaction_matrix, x)},
                {'type': 'ineq', 'fun': lambda x: x}]
        res = optimize.minimize(lsqmin, stoichiometric_coefficients_vector,
                                args=(reaction_rates_matrix, production_rates_vector), method='SLSQP',  # 'COBYLA',
                                constraints=cons, tol=1e-5,
                                options={'maxiter': 5000})

        # Storing the result correctly
        stoichiometric_coefficients_dict = {}
        counter = 0
        for reaction in self.products_by_reaction:
            stoichiometric_coefficients_dict[reaction] = {}
            for product in self.products_by_reaction[reaction]:
                if res.x[counter] > 1e-5:
                    stoichiometric_coefficients_dict[reaction][product] = round(res.x[counter], 5)
                counter += 1

        # TODO a supplementary step ensuring the elemental balance must be added
        # Checking the atoms conservation
        redo_balance = {}
        for index_reaction, reaction in enumerate(self.products_by_reaction):
            for index_element, element in enumerate(elements_per_reaction[reaction]):
                element_delta = total_elements_per_reaction[reaction][element] \
                                - sum(
                        [stoichiometric_coefficients_dict[reaction][spec] * self.composition_dict[spec][element] for
                         spec in
                         stoichiometric_coefficients_dict[reaction]])
                if abs(element_delta) > 1e-3:
                    if reaction not in redo_balance:
                        redo_balance[reaction] = {element: element_delta}
                    else:
                        redo_balance[reaction][element] = element_delta

        self.stoichiometric_coefficients_dict = stoichiometric_coefficients_dict

    def create_reactions(self):
        """Uses all the stored data to create Cantera reaction objects


        """

        import cantera as ct

        decomposition_reactions_list = []

        # Very ugly bugfix preventing the reactions of changing order in the mechanism
        # TODO must be made better !
        sorted_reactions_list = []
        if 'fuel' in self.Arrhenius_parameters_dict:
            sorted_reactions_list.append('fuel')
        self.precursors.sort()
        for spec in self.precursors:
            if spec in self.Arrhenius_parameters_dict:
                sorted_reactions_list.append(spec)

        # Create reaction
        for reaction in sorted_reactions_list:
            reactants_dict = {key: 1.0 for key in self.reactants_by_reaction[reaction]}
            products_dict = self.stoichiometric_coefficients_dict[reaction]
            # New elementary reaction
            new_reaction_object = ct.ElementaryReaction(reactants_dict, products_dict)

            if self.temperature_correction_Arrhenius:
                new_reaction_object.rate = ct.Arrhenius(self.Arrhenius_parameters_dict[reaction][0],
                                                        self.Arrhenius_parameters_dict[reaction][1],
                                                        self.Arrhenius_parameters_dict[reaction][2])
            else:
                new_reaction_object.rate = ct.Arrhenius(self.Arrhenius_parameters_dict[reaction][0],
                                                        0.0,
                                                        self.Arrhenius_parameters_dict[reaction][1])

            new_reaction_object.reversible = False
            decomposition_reactions_list.append(new_reaction_object)

        self.lumped_reactions_objects += decomposition_reactions_list

    def build_mechanism(self, name=None, discard_species=True):
        """Builds the Mechanism object

        Parameters
        ----------
        name :
            name of the new class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
        discard_species :
            if True, the species involved in pyrolysis are discarded,
            else only the pyrolysis reactions are (Default value = True)

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object after lumping

        """

        # Saving the fuel species
        if self.method == 'Pyrolysis' and self.fuel_species == ['Surrogate']:
            import ARCANE.lumping as lumping
            fuel_species_object = lumping.lump_non_isomer_species(self.mechanism, self.fuel_species_dict,
                                                                  lumped_species_name='Surrogate')
        else:
            fuel_species_object = [self.mechanism.species[self.mechanism.species_names.index(spec)]
                                   for spec in self.fuel_species_dict]

        self.mechanism.inert = ['N2']

        # The first step is to remove all the species flagged as decomposition intermediates
        if discard_species:
            trimmed_mechanism = self.mechanism.remove_species(self.reactants, auto_remove=True)
        # Or by default, remove all the decomposition reactions
        else:
            trimmed_mechanism = self.mechanism.remove_reactions(self.reactions, auto_remove=True)

        # Add the fuel_species once all its reactions have been removed
        trimmed_mechanism = trimmed_mechanism.add_species(fuel_species_object)

        new_mechanism = trimmed_mechanism.add_reaction(self.lumped_reactions_objects, name=name)

        return new_mechanism

    def reduce(self, name=None, method='Pyrolysis'):
        """Complete reduction of the pathways

        Parameters
        ----------
        name : str
            name of the new class :func:`~ARCANE.mechanisms.Mechanism` object (Default value = None)
        method : str
            Lumping method:
            - 'Pyrolysis': lumps the fuel species decomposition into HyChem-like reactions
            if the fuel contains several species, they will be lumped into an equivalent Surrogate species
            - 'Pyrolysis multi': lumps each fuel species decomposition an merge the resulting mechanisms
            - 'PAH': to come

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object after lumping

        """
        self.method = method
        if self.method == 'Pyrolysis' and len(self.fuel_species) > 1:
            self.fuel_species = ['Surrogate']
        else:
            self.fuel_species = list(self.fuel_species_dict.keys())

        global_fuel_species_dict = self.fuel_species_dict.copy()
        global_fuel_species = self.fuel_species.copy()

        self.lumped_reactions_objects = []

        for index, fuel_species in enumerate(global_fuel_species):

            if self.method == 'Pyrolysis multi':
                # Changing the base cases list with single fuel
                for case in self.cases_list:
                    case.change_case_input_parameter('fuel', f'X/{fuel_species}/1')

                # Run simulations for super root mechanism
                cases.run_cases(self.cases_list, self.mechanism)

                logger.info(f"\nPeforming the reduction for species {fuel_species}")

                self.fuel_species_dict = {fuel_species: 1}

            self.fuel_species = fuel_species

            # Loading the samples and relevant data
            self.load_samples()

            # Sorting the species and reactions according to their role in the mechanism
            self.get_pyrolysis_species()

            # Displaying relevant reactions
            # for index_dummy, reac in enumerate(self.reactions):
            #     print(sum(self.samples_database['net rates of progress'][index_dummy, :]))
            #     print(reac.equation)
            #quit()

            # Creating differentiated species source terms
            self.sorted_source_terms()

            # Computing the importance ratio of species involved in H abstraction
            self.compute_hydrogen_abstraction_importance()

            # Selecting the reactions
            sorted_list = sorted(zip(self.precursors_coefficients.values(),
                                     self.precursors_coefficients.keys()), key=lambda x: x[0])
            sorted_list.reverse()

            if self.use_threshold_for_reactions:
                # Selecting the species that are above the user input threshold
                self.selected_precursors = [spec for spec in self.precursors_coefficients
                                            if self.precursors_coefficients[spec]
                                            > self.reactants_importance_threshold[index]]
            else:
                if self.number_of_reactions[index] > 1:
                    self.selected_precursors = [value[1] for value in sorted_list[:self.number_of_reactions[index] - 1]]
                    self.reactants_importance_threshold[index] = self.precursors_coefficients[
                                                              self.selected_precursors[-1]] + 1e-60
                else:
                    self.selected_precursors = []
                    self.reactants_importance_threshold[index] = 100

            logger.info(f"\nThe possible H abstractions precursors and their coefficients are:")
            length_string = len(max([str(s[1]) for s in sorted_list], key=len))
            line_plotted = False
            for spec_value in sorted_list:
                if spec_value[0] < self.reactants_importance_threshold[index] and not line_plotted:
                    logger.info('-' * length_string + '----------------')
                    line_plotted = True
                logger.info(f"{spec_value[1]:{length_string}} => {spec_value[0]:.6E}")

            # Computing the rate coefficients of the lumped reactions
            self.compute_hydrogen_abstraction_rates()

            # If the sum of the selected precursors coefficients is greater than 1, the fuel homolysis is removed
            selected_precursors_coefficients = [self.precursors_coefficients[spec] for spec in self.selected_precursors]
            sum_coeffcients = sum(selected_precursors_coefficients)

            if sum_coeffcients < 1:
                self.compute_fuel_decomposition_rate()
            else:
                self.fuel_homolysis = False

            # Creating a single dictionary (one for each fuel component later) for the rate coefficients
            if self.fuel_homolysis:
                self.rate_coefficients = {'fuel': self.fuel_decomposition_rate_coefficient}
            else:
                self.rate_coefficients = {}

            for spec in self.selected_precursors:
                self.rate_coefficients[spec] = self.hydrogen_abstraction_rate_coefficients[spec]

            # Fitting the Arrhenius constants
            self.fit_rate_coeffcients()

            # Identifying the potential products
            self.identify_products()

            sorted_list = sorted(zip(self.products_coefficients.values(),
                                     self.products_coefficients.keys()), key=lambda x: x[0])
            sorted_list.reverse()

            # Selecting the products
            if self.use_threshold_for_products:
                self.selected_products = [spec for spec in self.products_coefficients
                                          if self.products_coefficients[spec] > self.products_importance_threshold[index]]
            else:
                self.selected_products = [value[1] for value in sorted_list[:self.number_of_products[index]]]
                self.products_importance_threshold[index] = self.products_coefficients[
                                                          self.selected_products[-1]] + 1e-60

            logger.info(f"\nThe possible products and their coefficients are:")
            length_string = len(max([str(s[1]) for s in sorted_list], key=len))
            line_plotted = False
            for spec_value in sorted_list:
                if spec_value[0] < self.products_importance_threshold[index] and not line_plotted:
                    logger.info('-' * length_string + '----------------')
                    line_plotted = True
                logger.info(f"{spec_value[1]:{length_string}} => {spec_value[0]:.6E}")

            if not self.selected_products:
                logger.error(
                    "ERROR: The products importance threshold is too high and therfore no species was selected.")
                quit()

            # Associating the products to each reaction
            composition_first_product = {}

            if self.fuel_homolysis:
                self.products_by_reaction = {'fuel': self.selected_products}
                self.reactants_by_reaction = {'fuel': [self.fuel_species]}
            else:
                self.products_by_reaction = {}
                self.reactants_by_reaction = {}

            # Adding the RH species for H-abstraction reactions
            revised_precursors_list = []
            for spec in self.selected_precursors:
                composition_first_product[spec] = self.species_objects(spec).composition
                if 'H' in composition_first_product[spec].keys():
                    composition_first_product[spec]['H'] += 1
                else:
                    composition_first_product[spec]['H'] = 1

                species_found = False
                for prod in self.mechanism.species_names:
                    if self.species_objects(prod).composition == composition_first_product[spec]:
                        self.products_by_reaction[spec] = [prod]
                        species_found = True

                if species_found:
                    self.products_by_reaction[spec] += self.selected_products
                    self.reactants_by_reaction[spec] = [self.fuel_species, spec]

                    # Removing duplicates
                    self.products_by_reaction[spec] = list(set(self.products_by_reaction[spec]))

                    revised_precursors_list.append(spec)
                else:
                    logger.info(
                            f"\nSpecies {spec} has been removed from the precursors because the corresponding RH cannot be found.")
                    del self.Arrhenius_parameters_dict[spec]

            # Checking if a product should be discarded because it creates elemental unbalance
            for reac in self.reactants_by_reaction:
                element_not_in_reaction = ['C', 'H', 'O']
                for spec in self.reactants_by_reaction[reac]:
                    for element in ['C', 'H', 'O']:
                        if self.composition_dict[spec][element] > 0 and element in element_not_in_reaction:
                            element_not_in_reaction.remove(element)

                # Removing
                revised_products_list = self.products_by_reaction[reac].copy()
                if element_not_in_reaction:
                    for element in element_not_in_reaction:
                        for prod in self.products_by_reaction[reac]:
                            if self.composition_dict[prod][element] != 0:
                                revised_products_list.remove(prod)

                self.products_by_reaction[reac] = revised_products_list
            # print(self.products_by_reaction)

            self.selected_precursors = revised_precursors_list

            # Solving the constrained linear system to find the stoichiometric coefficients
            self.solve_products_coefficients()

            logger.info(f"\n")

            # Creating the reaction objects
            self.create_reactions()

        # Displaying the results of the lumping
        logger.info(f"\nThe final lumped reactions are:")
        for reaction in self.lumped_reactions_objects:
            logger.info(reaction.equation)
            logger.info(reaction.rate)
        logger.info("\n")

        # Resetting the fuel species variables
        self.fuel_species_dict = global_fuel_species_dict
        self.fuel_species = global_fuel_species

        # Building the mechanism object
        new_mechanism = self.build_mechanism(name=name, discard_species=self.discard_species_in_build)

        head, tail = display.head_tail_charac(text=f'Mechanism {new_mechanism.name} created', size='small', style='@')
        logger.goodnews(head)

        # Changing the fuel back
        if self.method == 'Pyrolysis multi':
            # Changing the base cases list with single fuel
            for case in self.cases_list:
                case.change_case_input_parameter('fuel', self.fuel_species_dict)

        return new_mechanism
