# Module for generating a pyrolysis mechanism
import matplotlib.pyplot as plt

import cantera as ct

import ARCANE.cases as cases
import ARCANE.display as display
import ARCANE.sampling as sampling
import ARCANE.postproc as postproc
import ARCANE.tools as tools
import ARCANE.display as display
import ARCANE.analysis as analysis

import sys
import scipy.optimize as optimize
import numpy as np
import warnings
import scipy.optimize as opt

np.set_printoptions(threshold=sys.maxsize)
warnings.filterwarnings('ignore')

logger = display.Logger()
logger.set_log('logGrowth')
logger.create_log('logGrowth')


class Growth:

    def __init__(self, cases_list, mechanism, source, target, method='All', ranking='srate', condition='C>=6'):
        """
        """

        # Initialisation of important parameters
        self.method = method
        self.ranking = ranking
        self.mechanism = mechanism

        # Separating the cases for reduction from the complete list
        self.all_cases_list = cases_list
        self.cases_list = [case for case in cases_list]

        # Running the cases
        cases.run_cases(cases_list, mechanism)

        # Initial number of reactions
        self.nr = self.mechanism.nr

        # Number of samples for each case
        self.n_samples = 200000

        # Sets if the Arrhenius fit is done with the temperature correction
        self.temperature_correction_Arrhenius = False

        self.discard_species_in_build = True

        # List of final lumped reactions
        self.lumped_reactions_objects = []

        # Retrieve and store growth source and target species
        self.growth_source = [source]
        self.growth_target = [target]

        # Graph research thershold
        self.condition = condition

        # Truc
        self.source_species = source
        self.target_species = target

    def load_samples(self):
        """
        Loading the samples database

        :return:
        """

        self.samples_database = sampling.create_samples_database(self.cases_list, self.mechanism,
                                                                 n_samples=self.n_samples, time_normalised=True,
                                                                 integrated=False)

        self.samples_database = postproc.extend_data_dict(self.samples_database)

        self.n_samples = len(self.samples_database['grid'])

        return


    def get_all_growth_species_reactions(self):
        """
        Recover all species and reactions directly involved in any growth path identified through reaction graph analysis and roughly sorts species
        - growth_reactants : intermediates in the pathway. Will be removed in the lumped reactions
        - growth_species : precursors (reactants) or by-products in the final lumped reactions

        """

        paths, paths_weights, paths_reactions, paths_reactions_indices = analysis.growth_path(self.cases_list, self.mechanism, self.growth_source[0], self.growth_target[0], spec2plot=self.condition)

        # Recover species & reactions involved in any path
        growth_reactants = {s for s_list in paths for s in s_list}

        growth_reactions_indices = list({idx for sublist in paths_reactions_indices for idx in sublist})
        growth_reactions = [self.mechanism.reactions[idx] for idx in growth_reactions_indices]

        growth_species = []
        for reaction in growth_reactions:
            for reactant in reaction.reactants:
                growth_species.append(reactant)
            for product in reaction.products:
                growth_species.append(product)

        
        self.growth_species = list(set(growth_species) - growth_reactants)
        self.growth_reactants = list(growth_reactants - set(self.growth_source) - set(self.growth_target))
        self.growth_reactions = growth_reactions
        self.growth_reactions_indices = growth_reactions_indices

        return

    def get_main_growth_species_reactions(self):
        """
        Recover all species and reactions directly involved in any growth path identified through reaction graph analysis and roughly sorts species
        - growth_reactants : intermediates in the pathway. Will be removed in the lumped reactions
        - growth_species : precursors (reactants) or by-products in the final lumped reactions

        """

        paths, paths_weights, paths_reactions, paths_reactions_indices = analysis.growth_path(self.cases_list, self.mechanism, self.growth_source[0], self.growth_target[0], spec2plot=self.condition)

        # Recover species & reactions involved in most probable path
        main_idx = paths_weights.index(max(paths_weights))
        main_path = paths[main_idx]
        growth_reactants = {s for s in main_path}

        growth_reactions_indices = list({idx for idx in paths_reactions_indices[main_idx]})
        growth_reactions = [self.mechanism.reactions[idx] for idx in growth_reactions_indices]

        growth_species = []
        for reaction in growth_reactions:
            for reactant in reaction.reactants:
                growth_species.append(reactant)
            for product in reaction.products:
                growth_species.append(product)
        
        self.growth_species = list(set(growth_species) - growth_reactants)
        self.growth_reactants = list(growth_reactants - set(self.growth_source) - set(self.growth_target))
        self.growth_reactions = growth_reactions
        self.growth_reactions_indices = growth_reactions_indices

        return

    def compute_growth_reaction_importance(self):
        """
        """

        reactions_rates = self.samples_database['net rates of progress']
        grid = self.samples_database['grid']
        growth_reactions_rates = {}
        for reaction_index, reaction in enumerate(self.growth_reactions):
            global_reaction_index = self.growth_reactions_indices[reaction_index]
            reaction_rates = reactions_rates[global_reaction_index]
            reaction_rate = np.trapz(reaction_rates, grid)
            growth_reactions_rates[reaction] = reaction_rate
        
        self.growth_reactions_rates = growth_reactions_rates

        return

    def compute_growth_species_importance(self):
        '''
        '''
        molecular_weights = self.mechanism.ctmech.molecular_weights
        species_consumption_importance = dict.fromkeys(self.growth_species)
        species_production_importance = dict.fromkeys(self.growth_species)


        for species in self.growth_species:
            species_consumption_importance[species] = 0
            species_production_importance[species] = 0

        for reaction in self.growth_reactions_rates:
            rrate = self.growth_reactions_rates[reaction]
            reactants = reaction.reactants
            products = reaction.products
            for species in self.growth_species:
                if species in reactants:
                    if rrate > 0:
                        species_consumption_importance[species] += rrate
                    elif rrate < 0:
                        species_production_importance[species] -= rrate
                if species in products:
                    if rrate < 0:
                        species_consumption_importance[species] -= rrate
                    elif rrate > 0:
                        species_production_importance[species] += rrate
        
        grid = self.samples_database['grid']
        target_production_rates = self.samples_database['net production rates'][self.mechanism.ctmech.species_index(self.growth_source[0])] \
                                  * molecular_weights [self.mechanism.ctmech.species_index(self.growth_source[0])]
        target_production_rate = np.trapz(target_production_rates, grid)
        for species in self.growth_species:
            species_consumption_importance[species] *= molecular_weights [self.mechanism.ctmech.species_index(species)]
            species_consumption_importance[species] /= target_production_rate
            species_production_importance[species] *= molecular_weights [self.mechanism.ctmech.species_index(species)]
            species_production_importance[species] /= target_production_rate
        
        self.species_consumption_coefficients = species_consumption_importance
        self.species_production_coefficients = species_production_importance
        return



    def sorted_source_terms(self):
        """
        Sorts species by :
        - byproduct : not the target but possibly produced by lumped reaction (overally formed)
        - precursor : not the source, but possible reactant of the lumped reaction (overally consumed)
        """

        # Species and reactions of interest
        species = self.growth_source + self.growth_species + self.growth_target
        reactions = self.growth_reactions

        # Net reaction rates
        reaction_rates = self.samples_database['net rates of progress']

        # Sampling grid
        grid = self.samples_database['grid']


        growth_prod = {}
        growth_con = {}

        for spec in species:
            growth_prod[spec] = np.zeros([self.n_samples])
            growth_con[spec] = np.zeros([self.n_samples])

        # Loops over reactions regardless of its direction : 
        # if it mainly consumes the species of interest, it goes to consumption 
        # and vice-versa
        for reac_idx, reac in enumerate(reactions):
            global_reac_idx = self.growth_reactions_indices[reac_idx]
            reactants = list(reac.reactants.keys())
            products = list(reac.products.keys())
            # for sample in range(self.n_samples):
            #     rrate = reaction_rates[global_reac_idx, sample]
            #     for spec in species:
            #         global_spec_idx = self.mechanism.ctmech.species_index(spec)
            #         if spec in reactants:
            #             nu = ct.Kinetics.reactant_stoich_coeffs(self.mechanism.ctmech)[global_spec_idx, global_reac_idx]
            #             if rrate > 0 :
            #                 growth_con[spec][sample] += rrate * nu
            #             else:
            #                 growth_prod[spec][sample] -= rrate * nu
            #         elif spec in products:
            #             nu = ct.Kinetics.product_stoich_coeffs(self.mechanism.ctmech)[global_spec_idx, global_reac_idx]
            #             if rrate > 0:
            #                 growth_prod[spec][sample] += rrate * nu
            #             else:
            #                 growth_con[spec][sample] -= rrate * nu
            rrate = reaction_rates[global_reac_idx, :]
            sum_rrate = np.trapz(rrate, grid)
            for spec in species:
                global_spec_idx = self.mechanism.ctmech.species_index(spec)
                if spec in reactants:
                    nu = ct.Kinetics.reactant_stoich_coeffs(self.mechanism.ctmech)[global_spec_idx, global_reac_idx]
                    if sum_rrate > 0 :
                        growth_con[spec] += rrate * nu
                    else:
                        growth_prod[spec] -= rrate * nu
                elif spec in products:
                    nu = ct.Kinetics.product_stoich_coeffs(self.mechanism.ctmech)[global_spec_idx, global_reac_idx]
                    if sum_rrate > 0:
                        growth_prod[spec] += rrate * nu
                    else:
                        growth_con[spec] -= rrate * nu
                    
        self.growth_consumption_rates_dict = growth_con
        self.growth_production_rates_dict = growth_prod 

        return

    def compute_consumption_importance(self):
        alpha_dict = {}
        grid = self.samples_database['grid']

        # source_consumption_sum = 0
        # for source in self.growth_source:
        #     source_consumption = np.array(self.growth_consumption_rates_dict[source])
        #     source_consumption_sum += np.trapz(source_consumption, grid) * self.mechanism.ctmech.molecular_weights[self.mechanism.ctmech.species_index(source)]

        target_production_sum = 0
        for target in self.growth_target:
            target_production = np.array(self.growth_production_rates_dict[target])
            target_production_sum += np.trapz(target_production, grid) * self.mechanism.ctmech.molecular_weights[self.mechanism.ctmech.species_index(target)]

        for spec in self.growth_species:
            spec_consumption = np.array(self.growth_consumption_rates_dict[spec])
            spec_consumption_sum = np.trapz(spec_consumption, grid) * self.mechanism.ctmech.molecular_weights[self.mechanism.ctmech.species_index(spec)]
            if target_production_sum > 0:
                alpha = spec_consumption_sum / target_production_sum
            else:
                alpha = 0
            
            alpha_dict[spec] = alpha

        self.species_consumption_coefficients = alpha_dict

    def compute_production_importance(self):
        alpha_dict = {}
        grid = self.samples_database['grid']

        # source_consumption_sum = 0
        # for source in self.growth_source:
        #     source_consumption = np.array(self.growth_consumption_rates_dict[source])
        #     source_consumption_sum += np.trapz(source_consumption, grid) * self.mechanism.ctmech.molecular_weights[self.mechanism.ctmech.species_index(source)]

        target_production_sum = 0
        for target in self.growth_target:
            target_production = np.array(self.growth_production_rates_dict[target])
            target_production_sum += np.trapz(target_production, grid) * self.mechanism.ctmech.molecular_weights[self.mechanism.ctmech.species_index(target)]

        for spec in self.growth_species:
            spec_production = np.array(self.growth_production_rates_dict[spec])
            spec_production_sum = np.trapz(spec_production, grid) * self.mechanism.ctmech.molecular_weights[self.mechanism.ctmech.species_index(spec)]
            if target_production_sum > 0:
                alpha = spec_production_sum / target_production_sum
            else:
                alpha = 0

            alpha_dict[spec] = alpha


        self.species_production_coefficients = alpha_dict

    def fit_rate_coefficients(self):
        """

        :return:
        """
        
        Arrhenius_parameters_dict = dict.fromkeys(self.rate_coefficients.keys(), [])
        for spec in Arrhenius_parameters_dict:
            rate_coefficients_untreated = self.rate_coefficients[spec]
            # Pre-treat arrays to discard 0 values
            

            temperatures=np.array(self.samples_database[f"T"])
            # Mask for low concentration of pyrene
            too_low_concentration_mask = self.samples_database[f"c_{self.target_species}"] >  max(self.samples_database[f"c_{self.target_species}"]) / 10
            temperatures = temperatures[too_low_concentration_mask]
            

            temperatures[temperatures <= 0] = 'nan'
            rate_coeffs=rate_coefficients_untreated
            rate_coeffs = rate_coeffs[too_low_concentration_mask]  
            rate_coeffs[rate_coeffs <= 0] = 'nan'


            # Pre-treat arrays to discard NaN values
            temperatures = temperatures[~np.isnan(rate_coeffs)]
            rate_coeffs = rate_coeffs[~np.isnan(rate_coeffs)]
            # Pre-treat arrays to discard infs values
            temperatures = temperatures[~np.isinf(temperatures)]
            rate_coeffs = rate_coeffs[~np.isinf(rate_coeffs)]
                    
            #Checking 
            #print('rate_coeff', rate_coeffs)
            #print('temperatures', temperatures)
            self.temperature_correction_Arrhenius =  False
            # Fitting values
            if self.temperature_correction_Arrhenius:
                reverse_Arrh, pcov = opt.curve_fit(tools.calculate_lnk, temperatures, np.log(rate_coeffs))
            else:
                reverse_Arrh, pcov = opt.curve_fit(tools.calculate_lnk_simple, temperatures, np.log(rate_coeffs))


            # plt.plot(temperatures)
            # plt.show()
            plt.plot(1000/temperatures, np.log(rate_coeffs), 'o', label=spec)
            #plt.show()
            #plt.plot(1000 / temperatures, tools.calculate_lnk_simple(temperatures, *reverse_Arrh), label=spec + ' fit')
            #plt.show()
            #plt.plot(1000 / temperatures, tools.calculate_lnk(temperatures, *reverse_Arrh), label=spec + ' fit')
            plt.legend(loc="best", title="Data", frameon=False)
            #plt.show()
            print('revers Arhh', reverse_Arrh)
            #quit()
            # Retrieving real values
            reverse_Arrh[0] = np.exp(reverse_Arrh[0])
            if self.temperature_correction_Arrhenius:
                reverse_Arrh[2] = reverse_Arrh[2] * 8314.4621
            else:
                reverse_Arrh[1] = reverse_Arrh[1] * 8314.4621

            Arrhenius_parameters_dict[spec] = reverse_Arrh

        self.Arrhenius_parameters_dict = Arrhenius_parameters_dict
        print(self.Arrhenius_parameters_dict)
        # plt.legend()
        # plt.show()
        # quit()


    def create_reactions(self, append_reactions=True):
        """

        :return:
        """

        import cantera as ct

        reaction = "C2H2"
        decomposition_reactions_list = []
        # Create reaction
        #for reaction in self.Arrhenius_parameters_dict:
        reactants_dict = self.stoichiometric_coefficients_dict["1"][0]
        products_dict = self.stoichiometric_coefficients_dict["1"][1]
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
        #new_reaction_object.orders = {"C6H6":1, "C2H2":1, "H":1}
        new_reaction_object.orders = {"C10H7":1}

        self.lumped_reactions_objects += decomposition_reactions_list

    def build_mechanism(self, name=None, discard_species=False):
        """
        Builds the Mechanism object

        :param name:
        :param discard_species:
        :return:
        """

        # # Saving the fuel species
        # if self.method == 'Pyrolysis'and self.fuel_species == ['Surrogate']:
        #     import ARCANE.lumping as lumping
        #     fuel_species_object = lumping.lump_non_isomer_species(self.mechanism, self.fuel_species_dict,
        #                                                           lumped_species_name='Surrogate')
        # else:
        #     fuel_species_object = [self.mechanism.species[self.mechanism.species_names.index(spec)]
        #                            for spec in self.fuel_species_dict]

        self.mechanism.inert = ['N2']

        # The first step is to remove all the species flagged as decomposition intermediates
        #discard_species = False
        if discard_species:
            trimmed_mechanism = self.mechanism.remove_species(self.growth_reactants, auto_remove=True)
        # Or by default, remove all the decomposition reactions
        else:
            trimmed_mechanism = self.mechanism.remove_reactions(self.growth_reactions, auto_remove=True)
        print(self.growth_reactants)
        print(self.growth_reactions)
        # Add the fuel_species once all its reactions have been removed
       #trimmed_mechanism = trimmed_mechanism.add_species(fuel_species_object)

        new_mechanism = trimmed_mechanism.add_reaction(self.lumped_reactions_objects, name=name)
        #new_mechanism = trimmed_mechanism
        return new_mechanism

    def reduce(self, name =None):
        # Loading the samples and relevant data
        self.load_samples()
        # self.clip_samples()

        # Sorting the species and reactions according to their role in the mechanism
        if self.method=='All':
            self.get_all_growth_species_reactions()
        elif self.method=='Main':
            self.get_main_growth_species_reactions()
        else:
            raise NotImplementedError

        
        if self.ranking == "rrate":
            # Compute reactions rates of progress
            self.compute_growth_reaction_importance()
            # Compute species rates of formation / destruction through these reactions & normalize
            self.compute_growth_species_importance()
        elif self.ranking == "srate":
            # Compute production and consumption
            self.sorted_source_terms()
            # Normalize
            self.compute_production_importance()
            self.compute_consumption_importance()

        sorted_production_list = sorted(zip(self.species_production_coefficients.values(),
                                 self.species_production_coefficients.keys()), key=lambda x: x[0])

        sorted_consumption_list = sorted(zip(self.species_consumption_coefficients.values(),
                                 self.species_consumption_coefficients.keys()), key=lambda x: x[0]) 
            
        print(sorted_production_list)
        print(sorted_consumption_list)

        self.rate_coefficients = {}
        self.selected_precursors = self.species_consumption_coefficients.keys()
        print('selected precusros', self.selected_precursors)
        for spec in self.selected_precursors:
            #self.rate_coefficients[spec] = self.samples_database["c_C6H6"] #self.samples_database["cdot_C16H10"] / (  self.samples_database["c_C6H6"] )#* (self.samples_database["c_C2H2"])**1 )
            #self.rate_coefficients[spec] = self.samples_database["cdot_C16H10"] / (  self.samples_database["c_C6H6"] * (self.samples_database["c_C2H2"])**1 )#* (self.samples_database["c_C2H2"])**1 )
            self.rate_coefficients[spec] = self.samples_database[f"cdot_{self.target_species}"] / (  self.samples_database[f"c_{self.source_species}"])#* (self.samples_database[f"c_{spec}"])**1 )




        #print('rate coeffcient', self.rate_coefficients)
        mask = [True] * len(self.samples_database["c_C16H10"])#np.array(self.samples_database["c_C16H10"] > 1e-10)
        temp = np.array(self.samples_database['T'])[mask]
        conc_c16h10 = np.array(self.samples_database["c_C16H10"])[mask]
        conc_c6h6 = np.array(self.samples_database["c_C6H6"])[mask]
        conc_c2h2 = np.array(self.samples_database["c_C2H2"])[mask]

        w_c16h10 = np.array(self.samples_database["cdot_C16H10"])[mask]
        w_c6h6 = np.array(self.samples_database["cdot_C6H6"])[mask]
        w_c2h2 = np.array(self.samples_database["cdot_C2H2"])[mask]

        # plt.plot(temp,conc_c16h10 / max(conc_c16h10), label="C C16H10")
        # plt.plot(temp,conc_c6h6 / max(conc_c6h6) , label="C C6H6")
        # plt.plot(temp,conc_c2h2 /max(conc_c2h2), label="C C2H2")
        # plt.legend()
        # plt.show()
        
        # plt.plot(temp,w_c16h10 /max(w_c16h10) , label="w C16H10")
        # plt.plot(temp,w_c6h6 /max(w_c6h6), label="w C6H6")
        # plt.plot(temp,w_c2h2 / max(w_c2h2), label="W C2H2")
        # plt.legend()
        # # plt.plot(temp, np.log10(conc_c16h10))
        # # plt.plot(temp, np.log10(conc_c6h6))
        # # plt.plot(temp, np.log10(conc_c2h2))
        # plt.show()
        # quit()

        self.fit_rate_coefficients()

        # Impose stoichiometric coeffcients instead of using solve_product_coeff as in pyro
        # Smth like that :
        # 1 C6H6 + 5 C2H2 + H -> 1 C16H10 + 7 H
        # 1 C6H6 + 5 C2H2 + H -> 1 C16H10 + 7/2 H2

        stoichiometric_coefficients_dict = {}
        #stoichiometric_coefficients_dict = {"1":[{"C6H6":1, "C2H2":1, "H":1}, {"C16H10":1, "H":7 }]}
        stoichiometric_coefficients_dict = {"1":[{"C10H7":8}, {"C16H10":5, "H2":3}]}
        self.stoichiometric_coefficients_dict = stoichiometric_coefficients_dict

        # Creating the reaction 
        self.products_by_reaction = {}
        self.reactants_by_reaction = {}
        self.selected_products = [f"{self.target_species}", 'H']

        #for spec in self.selected_precursors:
        # self.products_by_reaction = {"C16H10":1, "H":7 }
        # self.reactants_by_reaction = {"C6H6":1, "C2H2":5, "H":1}

        self.products_by_reaction = {"C16H10":5, "H2":3 }
        self.reactants_by_reaction = {"C10H7":8}

        self.create_reactions() 

        # Creating the new mechanism
        new_mechanism = self.build_mechanism(name=name, discard_species=self.discard_species_in_build)
        new_mechanism = new_mechanism.remove_species(['INDENYL', 'C10H7CH2'])
        head, tail = display.head_tail_charac(text=f'Mechanism {name} created', size='small', style='@')


        logger.goodnews(head)
        
        return new_mechanism
