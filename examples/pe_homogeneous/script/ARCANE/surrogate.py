import numpy as np
import pandas as pd
from scipy.optimize import fsolve, minimize, LinearConstraint, curve_fit, OptimizeResult
import ARCANE.mechanisms as mechanisms
import ARCANE.kwdict as kwdict
import ARCANE.display as display
import ARCANE.chem_data as data
import ARCANE.AVBP as AVBP
import ARCANE.cases as cases

logger = display.Logger()
logger.set_log('logCompute')

kwdict = kwdict.Kwdict()


class Surrogate:

    def __init__(self, palette_mechanism, composition_dict, properties_dict=None, name=None, check_chem_db=False):
        """
        Generator for the Surrogate class

        Parameters
        ----------
        palette_mechanism : class :func:`ARCANE.mechanisms.Mechanism` object
            Mechanism object from ARCANE that is used for computation
        composition_dict : dict
            Dictionary containing the chromatography data or the species name and mole fraction
        properties_dict : dict
            Dictionary containing properties of the Surrogate object
        name : str
            Name of the Surrogate object
        check_chem_db : bool
            if true, species properties are calculated from ARCANE database (check chem_data.py)

        Created: 21/07/01 [QC]
        Last modified: 21/07/23 [TL]

        """

        # Initialising the attributes
        self.name = name
        self.mechanism = palette_mechanism
        self.composition = composition_dict
        self.check_chem_db = check_chem_db
        self.laminar_flame_velocity = pd.DataFrame({'temperature': [], 'pressure': [], 'phi': [], 'sl': []})
        self.ignition_delay_time = pd.DataFrame({'temperature': [], 'pressure': [], 'phi': [], 'tig': []})

        # Initialising the properties dict
        self.handled_properties = ['hydrogen_carbon_ratio', 'molecular_weight',
                                   'derived_cetane_number', 'threshold_sooting_index',
                                   'liquid_density', 'distillation_curve',
                                   'surface_tension', 'liquid_viscosity', 'low_heating_value']

        if properties_dict:
            self.properties = properties_dict
        else:
            self.properties = {}

        self.species_properties = {}

        # Defining the Hydrocarbon families to get data from chromatography (mole data)
        hydrocarbon_families = {"n-alkanes": {"formula": "C(n)H(2n+2)",
                                              "hydrogen number": lambda carbon_number: 2 * carbon_number + 2,
                                              "molecular weight": lambda
                                                  carbon_number: carbon_number * 12 + 2 * carbon_number + 2},
                                "i-alkanes": {"formula": "C(n)H(2n+2)",
                                              "hydrogen number": lambda carbon_number: 2 * carbon_number + 2,
                                              "molecular weight": lambda
                                                  carbon_number: carbon_number * 12 + 2 * carbon_number + 2},
                                "cyclo-alkanes": {"formula": "C(n)H(2n)",
                                                  "hydrogen number": lambda carbon_number: 2 * carbon_number,
                                                  "molecular weight": lambda
                                                      carbon_number: carbon_number * 12 + 2 * carbon_number},
                                "aromatics": {"formula": "C(n)H(2n-6)",
                                              "hydrogen number": lambda carbon_number: 2 * carbon_number - 6,
                                              "molecular weight": lambda
                                                  carbon_number: carbon_number * 12 + 2 * carbon_number - 6},
                                "di-aromatics": {"formula": "C(n)H(2n-12)",
                                                 "hydrogen number": lambda carbon_number: 2 * carbon_number - 12,
                                                 "molecular weight": lambda
                                                     carbon_number: carbon_number * 12 + 2 * carbon_number - 12},
                                }

        # Detecting the type of composition dictionary and checking it
        self.raw_composition = False
        #corrected_dict = {}
        for key in composition_dict:
            if key in hydrocarbon_families:
                self.raw_composition = True

            elif key in self.mechanism.species_names:
                if self.raw_composition:
                    logger.error("Both composition formulations cannot be used, either choose to express "
                                 "the dictionary in term of hydrocarbons families or as a surrogate")
                    quit()

            else:
                logger.error(f'Invalid entry {key} in the composition dict '
                             f'it should either be an hydrocarbon family or a species')
                quit()

        # Composition functions : raw_composition for chromatography
        if self.raw_composition:

            # Computing the quantities of interest for the detailed composition
            # determining the min and max carbon numbers

            families = list(composition_dict.keys())

            min_carbon = 100
            max_carbon = 0

            for key in composition_dict:
                carbon_numbers = list(composition_dict[key].keys())

                if carbon_numbers:
                    local_min = min(carbon_numbers)
                    if local_min < min_carbon:
                        min_carbon = min(carbon_numbers)

                    local_max = max(carbon_numbers)
                    if local_max > max_carbon:
                        max_carbon = max(carbon_numbers)

            carbons = list(range(min_carbon, max_carbon))

            # H/C ratio
            hydrogen_atoms = 0
            carbon_atoms = 0
            molecular_weight = 0

            for family in families:
                for carbon in carbons:
                    if carbon in composition_dict[family]:
                        carbon_atoms += carbon * composition_dict[family][carbon]

                        hydrogen_atoms += hydrocarbon_families[family]["hydrogen number"](carbon) \
                                          * composition_dict[family][carbon]

                        molecular_weight += hydrocarbon_families[family]["molecular weight"](carbon) \
                                            * carbon * composition_dict[family][carbon]

            self.properties['hydrogen_carbon_ratio'] = hydrogen_atoms / carbon_atoms
            self.properties['molecular_weight'] = molecular_weight / carbon_atoms

        # not raw_composition for data from species
        else:

            hydrogen_atoms = 0
            carbon_atoms = 0
            molecular_weight = 0

            for spec in composition_dict:
                # Check whether properties are known in the ARCANE chem database
                if check_chem_db:
                    self.species_properties[spec] = data.get_tpf_data(spec)

                # Calculation of H/C ratio from mechanism : mole averaging
                hydrogen_atoms += composition_dict[spec] * self.mechanism.ctmech.n_atoms(spec, 'H')
                carbon_atoms += composition_dict[spec] * self.mechanism.ctmech.n_atoms(spec, 'C')
                molecular_weight += composition_dict[spec] * self.mechanism.ctmech.molecular_weights[
                    self.mechanism.ctmech.species_index(spec)]

            self.properties['hydrogen_carbon_ratio'] = hydrogen_atoms / carbon_atoms
            self.properties['molecular_weight'] = molecular_weight

    def add_property(self, name, value):
        """
        Method to add a species in the property dictionary

        Parameters
        ----------

        name: str
            key name to add to the dictionary
        value:
            property value to be stored

        Created: 21/07/01 [TL]
        Last modified: 21/07/23 [TL]
        """
        if name in self.handled_properties:
            self.properties[name] = value
        else:
            logger.warning(f"The property {name} is not handled by the optim tool")
            logger.warning(f"Please give a property in the following list : {self.handled_properties}")

    

    def whole_range_composition_optimization(self, target_surrogate, list_target_properties, standard_deviation=np.zeros(0),
                                             weight=np.zeros(0), method='LMS', sl=False, tig=False, print_results=True):
        """

        Parameters
        ----------

        target_surrogate :

        list_target_properties :

        standard_deviation :

        weight :

        method :

        sl :

        tig :

        print_results :


        """

        # Store initial properties if results must be printed
        if print_results:
            list_old = []
            for prop in list_target_properties:
                if prop != 'distillation_curve':
                    self.calculate_property(prop)
                    list_old.append(self.properties[prop])

        ## Define initial conditions guesses on the whole range
        # Discrete mass fractions that can be taken
        nb_species = len(self.composition.keys())
        discrete_state = np.linspace(0, 1, nb_species + 1)
        initial_composition_list = []
        # Indices of the discrete state to be fixed for a given species
        indices = -np.ones(nb_species, dtype=int)

        # Loop on all the possible state
        for i in range((nb_species + 1) ** nb_species):
            # stored data for vector
            vector = np.zeros(nb_species)

            # update the species indices to cover all the possible state
            for j in range(nb_species):
                if j < nb_species - 1:
                    if (i % ((nb_species + 1) ** (nb_species - 1 - j))) == 0:
                        indices[j] += 1
                else:
                    indices[j] += 1
                if (indices[j] % (nb_species + 1)) == 0:
                    indices[j] = 0

                # for each species, update the state
                vector[j] = discrete_state[indices[j]]

            # if the sum of all the mass fraction is one, than the state is possible
            if np.sum(vector) == 1: initial_composition_list.append(vector)

        # Perform optimization for each initial composition
        result = {}
        minimum = 1
        optim_mole_fraction = np.zeros(nb_species)
        for index, composition in enumerate(initial_composition_list):
            self.composition_optimization(target_surrogate, list_target_properties, initial_guesses=composition,
                                     standard_deviation=standard_deviation, weight=weight, method=method, sl=sl,
                                     tig=tig, print_results=False,result=result)
            if index > 0 :
                if result['fun'] < minimum :
                    minimum = result['fun']
                    for ind, mole_fraction in enumerate(result['x']):
                        if mole_fraction < 0.01:
                            result['x'][ind] = 0
                    optim_mole_fraction[:] = result['x'][:] / np.sum(result['x'])
            else:
                minimum = result['fun']
                for ind, mole_fraction in enumerate(result['x']):
                    if mole_fraction < 0.01:
                        result['x'][ind] = 0
                optim_mole_fraction[:] = result['x'][:] / np.sum(result['x'])

        for index, spec in enumerate(self.composition):
            self.composition[spec] = optim_mole_fraction[index]

        # Print the results if asked
        if print_results:
            print(" ")
            print("#####################   Optimization Results   #######################")
            for i, prop in enumerate(list_target_properties):
                if prop != 'distillation_curve':
                    print(
                        f"{prop}: {list_old[i]} => {self.properties[prop]}, objective was {target_surrogate.properties[prop]}")
            print(f"Optimal composition (mole frac.) found is : {self.composition}")
            print(" ")

        return


    def composition_optimization(self, target_surrogate, list_target_properties, initial_guesses=np.zeros(0),
                                 standard_deviation=np.zeros(0), weight=np.zeros(0), method='LMS', sl=False, tig=False,
                                 print_results=True, result={}):
        """
        Composition optimization for a given surrogate (fixed species) for
        a list of target properties compared to a target surrogate

        Parameters
        ----------

        target_surrogate:
            Surrogate object to mimic
        list_target_properties:
            list of target properties to match
        initial_guesses:
            initial guesses for mole fraction, if not given all mole fractions are equal
        standard_deviation:
            np.array with the standard deviation expressed in % of the target mean value,
            0.2 by default
        weight:
            np.array for the weight added for each properties in the cost function, 1 for each fields
            by default
        method:
            formulation of the cost function, 'LMS' (least mean square) or 'EF' (exponential formulation)
        sl:
            if True, laminar flame speed computation is added in the cost function
        tig:
            if True, autoignition delay time is added in the cost function
        print_results:
            if True, results of the optimization are written in the terminal
        result:
            argument to get the result of the optimization back


        Created: 21/07/01 [TL]
        Last modified: 21/07/23 [TL]
        """

        # Internal definition of the cost function for optimization
        def objective_function(x, list_target_properties, varying_surrogate, target_surrogate, standard_deviation,
                               weight, method, sl, tig):
            f = 0
            for index, spec in enumerate(self.composition):
                self.composition[spec] = x[index]

            for j, prop in enumerate(list_target_properties):
                varying_surrogate.calculate_property(prop)
                if prop != 'distillation_curve':
                    if method == 'LMS':
                        f += weight[j] * ((varying_surrogate.properties[prop] - target_surrogate.properties[prop]) /
                                          target_surrogate.properties[prop]) ** 2
                    elif method == 'EF':
                        exp_argument = -(varying_surrogate.properties[prop] - target_surrogate.properties[prop]) ** 2 / \
                                       (2 * (standard_deviation[j] * target_surrogate.properties[prop]) ** 2)
                        f += weight[j] * (1 - np.exp(exp_argument)) ** 2
                    else:
                        logger.error(f"Objective formulation {method} is not known, choose LMS or EF")
                        quit()
                else:
                    varying_surrogate.match_abscissa(target_surrogate, prop)
                    temp_f = 0
                    for i in range(len(varying_surrogate.properties[prop][0])):
                        if method == 'LMS':
                            temp_f += ((varying_surrogate.properties[prop][1][i] - target_surrogate.properties[prop][1][
                                i]) / target_surrogate.properties[prop][1][i]) ** 2
                        elif method == 'EF':
                            exp_argument = - (varying_surrogate.properties[prop][1][i] -
                                              target_surrogate.properties[prop][1][i]) ** 2 / \
                                           (2 * (standard_deviation[j] * np.mean(
                                               target_surrogate.properties[prop][1])) ** 2)
                            temp_f += (1 - np.exp(exp_argument)) ** 2
                        else:
                            logger.error(f"Objective formulation {method} is not known, choose LMS or EF")
                            quit()

                    f += weight[j] / len(varying_surrogate.properties[prop][0]) * temp_f

            if sl:
                varying_surrogate.calc_laminar_flame_velocity_from_target(target_surrogate)
                temp_f = 0
                for i in range(len(target_surrogate.laminar_flame_velocity.values)):
                    if method == 'LMS':
                        temp_f += ((varying_surrogate.laminar_flame_velocity['sl'][i] -
                                    target_surrogate.laminar_flame_velocity['sl'][i])
                                   / target_surrogate.laminar_flame_velocity['sl'][i]) ** 2
                    elif method == 'EF':
                        exp_argument = - (varying_surrogate.laminar_flame_velocity['sl'][i] -
                                          target_surrogate.laminar_flame_velocity['sl'][i]) ** 2 / \
                                       (2 * (standard_deviation[len(list_target_properties)] * np.mean(
                                           target_surrogate.laminar_flame_velocity['sl'])) ** 2)
                        temp_f += (1 - np.exp(exp_argument)) ** 2
                    else:
                        logger.error(f"Objective formulation {method} is not known, choose LMS or EF")
                        quit()

                f += weight[len(list_target_properties)] / len(target_surrogate.laminar_flame_velocity.values) * temp_f

            if tig:
                varying_surrogate.calc_ignition_delay_time_from_target(target_surrogate)
                temp_f = 0
                for i in range(len(target_surrogate.ignition_delay_time.values)):
                    if method == 'LMS':
                        temp_f += ((varying_surrogate.ignition_delay_time['tig'][i] -
                                    target_surrogate.ignition_delay_time['tig'][i])
                                   / target_surrogate.ignition_delay_time['tig'][i]) ** 2
                    elif method == 'EF':
                        exp_argument = - (varying_surrogate.ignition_delay_time['tig'][i] -
                                          target_surrogate.ignition_delay_time['tig'][i]) ** 2 / \
                                       (2 * (standard_deviation[len(list_target_properties) + 1] * np.mean(
                                           target_surrogate.ignition_delay_time['tig'])) ** 2)
                        temp_f += (1 - np.exp(exp_argument)) ** 2
                    else:
                        logger.error(f"Objective formulation {method} is not known, choose LMS or EF")
                        quit()

                f += weight[len(list_target_properties)] / len(target_surrogate.ignition_delay_time.values) * temp_f
            return f

        # Check before calculation : composition is by species for the varying surrogate
        if self.raw_composition:
            logger.error('Composition should given by species to perform composition optimization')
            quit()

        # Check before calculation : target properties in input must be defined for target surrogate
        for prop in list_target_properties:
            if not prop in target_surrogate.properties:
                logger.error(f"The property {prop} is not defined for the target surrogate")
                quit()

        # Store initial properties if results must be printed
        if print_results:
            list_old = []
            for prop in list_target_properties:
                if prop != 'distillation_curve':
                    self.calculate_property(prop)
                    list_old.append(self.properties[prop])

        # Composition initialisation and constraints
        number_species = len(self.composition.keys())
        x_l_0 = np.zeros(number_species)
        if len(initial_guesses) == 0:
            x_l_0[:] = 1.0 / number_species
        else:
            if len(initial_guesses) == len(self.composition.keys()):
                for index, spec in enumerate(self.composition):
                    x_l_0[index] = initial_guesses[index]
            else:
                logger.error("The length of the initial guesses should match the number of species in the surrogate")
                quit()
        constraint = LinearConstraint(np.ones(number_species), lb=1.0, ub=1.0)
        bounds = [(0.0, 1.0) for i in range(number_species)]

        # Defined a default standard_deviation if user has not defined it (useful for EF formulation)
        if len(standard_deviation) == 0:
            standard_deviation = 0.2 * np.ones(len(list_target_properties) + 2)
        if len(weight) == 0:
            weight = np.ones(len(list_target_properties) + 2)

        # Call the optimization function, if tig or sl are computed the initial step gradient (eps) must be sufficient
        # to get a different computation
        if tig or sl:
            res = minimize(objective_function, x0=x_l_0,
                           args=(
                           list_target_properties, self, target_surrogate, standard_deviation, weight, method, sl, tig),
                           constraints=constraint, bounds=bounds, method='SLSQP', options={'eps': 0.015})
        else:
            res = minimize(objective_function, x0=x_l_0,
                           args=(
                           list_target_properties, self, target_surrogate, standard_deviation, weight, method, sl, tig),
                           constraints=constraint, bounds=bounds)

        # Get the result in a variable
        result['fun'] = res['fun']
        result['x'] = res['x']

        # Update the varying surrogate composition with the results of the optimization
        for index, mole_fraction in enumerate(res['x']):
            if mole_fraction < 0.01:
                res['x'][index] = 0
        res['x'][:] = res['x'][:] / np.sum(res['x'])

        for index, spec in enumerate(self.composition):
            self.composition[spec] = res['x'][index]

        # Print the results if asked
        if print_results:
            print(" ")
            print("#####################   Optimization Results   #######################")
            for i, prop in enumerate(list_target_properties):
                if prop != 'distillation_curve':
                    print(
                        f"{prop}: {list_old[i]} => {self.properties[prop]}, objective was {target_surrogate.properties[prop]}")
            print(f"Optimal composition found (mole frac.) is : {self.composition}")
            print(" ")

        return

    def match_abscissa(self, target_surrogate, prop):
        """
        Transform the data of the property to match the data of the target surrogate

        Parameters
        ----------

        target_surrogate:
            Surrogate object to match
        prop:
            name of the property

        :return: None

        Created: 21/07/01 [TL]
        Last modified: 21/07/23 [TL]
        """

        # Fourth order polinomial used to match the target property
        def fit_function(x, a, b, c, d, e):
            return a * x ** 4 + b * x ** 3 + c * x ** 2 + d * x + e

        # Check if the property exists in both surrogate objects
        if (not prop in self.properties) or (not prop in target_surrogate.properties):
            logger.warning(f"The property {prop} is not provided for one of the surrogate")
            return
        # Check if the property to match has the right format
        if not isinstance(self.properties[prop], list) or not isinstance(target_surrogate.properties[prop], list):
            logger.error(f"The property {prop} given in fit_function should have the following format [[x],[y]]")
            return
        if len(self.properties[prop]) != 2 or len(target_surrogate.properties[prop]) != 2:
            logger.error(f"The property {prop} given in fit_function should have the following format [[x],[y]]")
            return

        # Use scipy curve_fit to match a polinomial to the data
        popt, _ = curve_fit(fit_function, self.properties[prop][0],
                            self.properties[prop][1])

        # Get the matched coefficient for the polinomial
        a, b, c, d, e = popt

        # Update the field in the surrogate to match the target data points
        abscissa = np.array(target_surrogate.properties[prop][0])
        ordinate = fit_function(abscissa, a, b, c, d, e)
        del self.properties[prop]
        self.properties[prop] = [abscissa, ordinate]

        return

    def calculate_property(self, prop, x=np.zeros(0)):
        """
        Calculates a surrogate property with the appropriate mixing rule and updates composition

        Parameters
        ----------

        prop:
            Name of the property
        x:
            np.array for mole fractions. If given, update the surrogate composition with the mole fractions given

        :return: None

        Created: 21/07/01 [TL]
        Last modified: 21/07/23 [TL]
        """
        # Check before calculation : composition is by species
        if self.raw_composition:
            logger.error('Composition should given by species to perform composition optimization')
            quit()

        # Update surrogate composition
        if len(x) > 0:
            for index, spec in enumerate(self.composition):
                self.composition[spec] = x[index]

        # Calculates the property with the right mixing rule
        if prop == 'hydrogen_carbon_ratio':
            hydrogen_atoms = 0
            carbon_atoms = 0
            for spec in self.composition:
                # Calculation of H/C ratio from mechanisms, mole averaging
                hydrogen_atoms += self.composition[spec] * self.mechanism.ctmech.n_atoms(spec, 'H')
                carbon_atoms += self.composition[spec] * self.mechanism.ctmech.n_atoms(spec, 'C')

            self.properties['hydrogen_carbon_ratio'] = hydrogen_atoms / carbon_atoms
            return

        elif prop == 'molecular_weight':
            molecular_weight = 0
            for spec in self.composition:
                molecular_weight += self.composition[spec] * self.mechanism.ctmech.molecular_weights[
                    self.mechanism.ctmech.species_index(spec)]

            self.properties['molecular_weight'] = molecular_weight
            return

        elif prop == 'liquid_density':
            density = 0
            for spec in self.composition:
                density += self.composition[spec] * self.species_properties[spec]['liquid_density']

            self.properties['liquid_density'] = density
            return

        elif prop == 'threshold_sooting_index':
            tsi = 0
            for spec in self.composition:
                tsi += self.composition[spec] * self.species_properties[spec]['threshold_sooting_index']

            self.properties['threshold_sooting_index'] = tsi
            return

        elif prop == 'derived_cetane_number':
            dcn = 0
            for spec in self.composition:
                dcn += self.composition[spec] * self.species_properties[spec]['derived_cetane_number']

            self.properties['derived_cetane_number'] = dcn
            return

        elif prop == 'surface_tension':
            self.calc_surface_tension()
            return

        elif prop == 'distillation_curve':
            self.calc_distillation_curve()
            return

        elif prop == 'low_heating_value':
            self.calc_low_heating_value()
            return

        elif prop == 'liquid_viscosity':
            self.calc_liquid_viscosity()
            return

        else:
            logger.error(f"The property {prop} is not implemented in the function")
            quit()
            return

    def define_palette(self, list_of_palette_species):
        """
        Defines the species that can be used for the definition of the surrogate

        Parameters
        ----------

        list_of_palette_species:
            list of palette species to be used

        :return: None

        Created: 21/07/01 [QC]
        """

        # Checking if the species are in the mechanism
        valid_species_bool = [spec in self.mechanism for spec in list_of_palette_species]

        if not all(valid_species_bool):
            invalid_species = [spec for index, spec in enumerate(list_of_palette_species)
                               if not valid_species_bool[index]]
            valid_species = [spec for index, spec in enumerate(list_of_palette_species)
                             if valid_species_bool[index]]
            invalid_species_string = (', ').join(invalid_species)
            logger.error(f"{invalid_species_string} should be in the mechanism.")
            logger.error(f"Valid species have been included in the palette.")

        else:
            valid_species = list_of_palette_species

        self.palette_species = valid_species

        return


    def calc_low_heating_value(self):
        """
        Calculates the low heating value for the surrogate, adapted from :
        https://cantera.org/examples/jupyter/thermo/heating_value.ipynb.html

        Created: 21/07/01 [TL]
        """

        # Define the fuel composition in a Cantera-friendly form
        fuel = ""
        for spec in self.composition:
            fuel += f"{spec}:{format(self.composition[spec], '.4f')},"
        fuel = fuel[:-1]

        # Compute LHV
        self.mechanism.ctmech.TP = 298, 101325
        self.mechanism.ctmech.set_equivalence_ratio(1.0, fuel, 'O2:1.0')
        h1 = self.mechanism.ctmech.enthalpy_mass
        y_fuel = 0
        for spec in self.composition:
            y_fuel += self.mechanism.ctmech[spec].Y[0]

        # complete combustion products
        y_products = {'CO2': self.mechanism.ctmech.elemental_mole_fraction('C'),
                      'H2O': 0.5 * self.mechanism.ctmech.elemental_mole_fraction('H'),
                      'N2': 0.5 * self.mechanism.ctmech.elemental_mole_fraction('N')}

        self.mechanism.ctmech.TPX = None, None, y_products
        # Y_H2O = self.mechanism.ctmech['H2O'].Y[0]
        h2 = self.mechanism.ctmech.enthalpy_mass
        LHV = -(h2 - h1) / y_fuel
        # HHV = -(h2 - h1 + (h_liquid - h_gas) * Y_H2O) / Y_fuel

        self.properties['low_heating_value'] = LHV / 10 ** 6
        return LHV

    def calc_surface_tension(self, temperature=350):
        """
        Calculates the surface tension for a mixture, simplified MacLeod and Sugden

        Parameters
        ----------

        temperature: float
            temperature at which the surface tension is computed

        Created: 21/07/10 [TL]
        """
        sigma = 0
        for index, spec in enumerate(self.composition):
            # MacLeod and Sugden : sigma**0.25 = somme(Parachor_i * (xi*rho_mol_liq - yi*rho_mol_gas)
            # Parachor calculation (P = sigma**0.25/(rho_mol_liq - rho_mol_gas), rho_mol_gas neglected)
            # Then reduced to sigma**0.25 = somme(sigma_i**0.25 * xi)
            if not 'surface_tension' in self.species_properties[spec]:
                self.species_properties[spec]['surface_tension'] = AVBP.species_surface_tension(
                    self.species_properties[spec], calc_temperature=temperature)
            sigma += self.composition[spec] * (self.species_properties[spec]['surface_tension'] ** 0.25)

        self.properties['surface_tension'] = sigma ** 4
        return

    def calc_liquid_viscosity(self, temperature=253):
        """
        Calculates the liquid viscosity for a mixture, Gambill correlation

        Parameters
        ----------

        temperature:
            temperature at which the liquid viscosity is computed, -20°C by default

        Created: 21/07/10 [TL]
        """
        mu = 0
        for index, spec in enumerate(self.composition):
            if not 'liquid_viscosity' in self.species_properties[spec]:
                visc_list = AVBP.species_liquid_viscosity(self.species_properties[spec])
                self.species_properties[spec]['liquid_viscosity'] = visc_list[int(temperature / 10)]
            # Gambill 1959 correlation, to be improved (Chevron correl seems interesting but requires kinematic visco)
            mu += self.composition[spec] * (self.species_properties[spec]['liquid_viscosity'] ** (1 / 3))
        self.properties['liquid_viscosity'] = mu ** 3
        return

    def calc_laminar_flame_velocity(self, temperature, pressure, phi, oxidizer="X/O2/1.0/N2/3.76",
                                    method="approximate"):
        """
        Calculates the laminar flame velocity for a mixture, Ji correlation by default
        (https://doi.org/10.1016/j.proci.2010.06.085)

        Parameters
        ----------

        temperature:
            temperature string in ARCANE case definition style
        pressure:
            pressure string in ARCANE case definition style
        phi:
            phi string in ARCANE case definition style
        oxidizer:
            oxidizer string in ARCANE case definition style
        method:
            if "approximate", computation with Ji correlation. Otherwise, by Cantera

        Created: 21/07/15 [TL]
        """

        # Computation with Ji correlation
        if method == "approximate":

            # Local initialization
            temperature_adiab_mix = np.zeros(0)
            nx = np.zeros(0)  # product of number of products moles (C02+H20+N2) and mole fraction composition
            log_sli = np.zeros(0)
            sl = np.zeros(0)
            stored_cond_thermo = np.zeros(0)

            # Loop on composition to calculate each contribution to sl
            for i, spec in enumerate(self.composition):
                # Calculate each component laminar flame velocity
                fuel = f"X/{spec}/1.0"
                caselist = []
                caselist.extend(cases.create_case(reactor="C1DP",
                                                  mechanism=self.mechanism,
                                                  fuel=fuel,
                                                  oxidizer=oxidizer,
                                                  pressure=str(pressure),
                                                  temperature=str(temperature),
                                                  phi=str(phi)))
                cases.run_cases(caselist, self.mechanism, overwrite=False)

                # Initialize variables for each cases
                if i == 0:
                    temperature_adiab_mix = np.zeros(len(caselist))
                    nx = np.zeros(
                        len(caselist))  # product of number of products moles (C02+H20+N2) and mole fraction composition
                    log_sli = np.zeros(len(caselist))
                    sl = np.zeros(len(caselist))
                    stored_cond_thermo = np.zeros([3, len(caselist)])

                # For each case, get the data of interest
                for j, case in enumerate(caselist):
                    sl = case.extract_quantity(self.mechanism, 'Laminar flame speed')

                    # Calculate adiabatic temperature for each temperature
                    fuel = f"{spec}:1.0"
                    data_dict_from_case = case.data_dict(self.mechanism)
                    stored_cond_thermo[0, j] = data_dict_from_case['Temperature'][0]
                    stored_cond_thermo[1, j] = data_dict_from_case['Pressure'][0]
                    stored_cond_thermo[2, j] = case.parameters['phi']
                    self.mechanism.ctmech.TP = data_dict_from_case['Temperature'][0], data_dict_from_case['Pressure'][0]
                    self.mechanism.ctmech.set_equivalence_ratio(case.parameters['phi'], fuel, 'O2:1.0,N2:3.76')
                    number_products_moles = case.parameters['phi'] * self.mechanism.ctmech.n_atoms(spec, 'C') \
                                            + case.parameters['phi'] * self.mechanism.ctmech.n_atoms(spec, 'H') / \
                                            2 + 3.76 * (self.mechanism.ctmech.n_atoms(spec, 'C') \
                                                        + self.mechanism.ctmech.n_atoms(spec, 'H') / 4)
                    self.mechanism.ctmech.equilibrate('HP')
                    temperature_adiab = self.mechanism.ctmech.T

                    # Compute the contribution of each case
                    temperature_adiab_mix[j] += self.composition[spec] * number_products_moles * temperature_adiab
                    nx[j] += self.composition[spec] * number_products_moles
                    log_sli[j] += self.composition[spec] * number_products_moles * temperature_adiab * np.log(sl)

            #  Normalize the data
            temperature_adiab_mix = temperature_adiab_mix / nx
            sl = np.exp(log_sli / temperature_adiab_mix / nx)

            # Store the data computed
            for i in range(len(stored_cond_thermo[1, :])):
                add_row = pd.Series(data=[stored_cond_thermo[0, i],
                                          stored_cond_thermo[1, i],
                                          stored_cond_thermo[2, i],
                                          sl[i]],
                                    index=self.laminar_flame_velocity.columns,
                                    name=self.laminar_flame_velocity.shape[0])
                self.laminar_flame_velocity = self.laminar_flame_velocity.append(add_row)

        # Computation with Cantera
        else:
            # Fuel definition in a ARCANE friendly style
            fuel = "X/"
            for spec in self.composition:
                fuel += f"{spec}/{format(self.composition[spec], '.4f')}/"
            fuel = fuel[:-1]

            # Define cases to run
            caselist = []
            caselist.extend(cases.create_case(reactor="C1DP",
                                              mechanism=self.mechanism,
                                              fuel=fuel,
                                              oxidizer=oxidizer,
                                              pressure=str(pressure),
                                              temperature=str(temperature),
                                              phi=str(phi)))

            cases.run_cases(caselist, self.mechanism, overwrite=False)

            # Store the data
            for case in caselist:
                data_dict_from_case = case.data_dict(self.mechanism)
                add_row = pd.Series(data=[data_dict_from_case['Temperature'][0],
                                          data_dict_from_case['Pressure'][0],
                                          case.parameters['phi'],
                                          case.extract_quantity(self.mechanism, 'Laminar flame speed')],
                                    index=self.laminar_flame_velocity.columns,
                                    name=self.laminar_flame_velocity.shape[0])
                self.laminar_flame_velocity = self.laminar_flame_velocity.append(add_row)

        return

    def calc_laminar_flame_velocity_from_target(self, target, oxidizer="X/O2/1.0/N2/3.76", method="approximate"):
        """
        Calculates the laminar flame velocity for a mixture based on the thermo condition of the target,
        Ji correlation applied by default (https://doi.org/10.1016/j.proci.2010.06.085)

        Parameters
        ----------

        target:
            Surrogate object to match
        oxidizer:
            oxidizer string in ARCANE case definition style
        method:
            if "approximate", computation with Ji correlation. Otherwise, by Cantera

        Created: 21/07/15 [TL]
        """

        # Erase data in the dataframe
        self.laminar_flame_velocity.drop(self.laminar_flame_velocity.index, inplace=True)

        # Formatting thermo conditions to ARCANE cases definition
        added_pressure = []
        added_temperature = []
        added_phi = []
        pressure = ""
        temperature = ""
        phi = ""
        for i in range(len(target.laminar_flame_velocity['temperature'])):
            if not target.laminar_flame_velocity['temperature'][i] in added_temperature:
                added_temperature.append(target.laminar_flame_velocity['temperature'][i])
                temperature += f'{target.laminar_flame_velocity["temperature"][i]}-'
            if not target.laminar_flame_velocity['pressure'][i] in added_pressure:
                added_pressure.append(target.laminar_flame_velocity['pressure'][i])
                pressure += f'{target.laminar_flame_velocity["pressure"][i]}-'
            if not target.laminar_flame_velocity['phi'][i] in added_phi:
                added_phi.append(target.laminar_flame_velocity['phi'][i])
                phi += f'{target.laminar_flame_velocity["phi"][i]}-'
        phi = phi[:-1]
        temperature = temperature[:-1]
        pressure = pressure[:-1]

        # Compute laminar flame speeds with Ji correlation
        if method == "approximate":

            # Local initialization
            temperature_adiab_mix = np.zeros(0)
            nx = np.zeros(0)  # product of number of products moles (C02+H20+N2) and mole fraction composition
            log_sli = np.zeros(0)
            sl = np.zeros(0)
            stored_cond_thermo = np.zeros(0)

            # Loop on composition to calculate each contribution to sl
            for i, spec in enumerate(self.composition):

                # Calculates each component laminar flame velocity
                fuel = f"X/{spec}/1.0"
                caselist = []
                caselist.extend(cases.create_case(reactor="C1DP",
                                                  mechanism=self.mechanism,
                                                  fuel=fuel,
                                                  oxidizer=oxidizer,
                                                  pressure=str(pressure),
                                                  temperature=str(temperature),
                                                  phi=str(phi)))
                cases.run_cases(caselist, self.mechanism, overwrite=False)

                # Initialize variables for each cases
                if i == 0:
                    temperature_adiab_mix = np.zeros(len(caselist))
                    nx = np.zeros(
                        len(caselist))  # product of number of products moles (C02+H20+N2) and mole fraction composition
                    log_sli = np.zeros(len(caselist))
                    sl = np.zeros(len(caselist))
                    stored_cond_thermo = np.zeros([3, len(caselist)])


                for j, case in enumerate(caselist):
                    sl = case.extract_quantity(self.mechanism, 'Laminar flame speed')

                    # Calculate adiabatic temperature for each condition defined
                    fuel = f"{spec}:1.0"
                    data_dict_from_case = case.data_dict(self.mechanism)
                    stored_cond_thermo[0, j] = data_dict_from_case['Temperature'][0]
                    stored_cond_thermo[1, j] = data_dict_from_case['Pressure'][0]
                    stored_cond_thermo[2, j] = case.parameters['phi']
                    self.mechanism.ctmech.TP = data_dict_from_case['Temperature'][0], data_dict_from_case['Pressure'][0]
                    self.mechanism.ctmech.set_equivalence_ratio(case.parameters['phi'], fuel, 'O2:1.0,N2:3.76')
                    number_products_moles = case.parameters['phi'] * self.mechanism.ctmech.n_atoms(spec, 'C') \
                                            + case.parameters['phi'] * self.mechanism.ctmech.n_atoms(spec, 'H') / \
                                            2 + 3.76 * (self.mechanism.ctmech.n_atoms(spec, 'C') \
                                                        + self.mechanism.ctmech.n_atoms(spec, 'H') / 4)
                    self.mechanism.ctmech.equilibrate('HP')
                    temperature_adiab = self.mechanism.ctmech.T

                    # Compute the contribution of each case
                    temperature_adiab_mix[j] += self.composition[spec] * number_products_moles * temperature_adiab
                    nx[j] += self.composition[spec] * number_products_moles
                    log_sli[j] += self.composition[spec] * number_products_moles * temperature_adiab * np.log(sl)

            # Normalize the data
            temperature_adiab_mix = temperature_adiab_mix / nx
            sl = np.exp(log_sli / temperature_adiab_mix / nx)

            # Store the data
            for i in range(len(stored_cond_thermo[1, :])):
                add_row = pd.Series(data=[stored_cond_thermo[0, i],
                                          stored_cond_thermo[1, i],
                                          stored_cond_thermo[2, i],
                                          sl[i]],
                                    index=self.laminar_flame_velocity.columns,
                                    name=self.laminar_flame_velocity.shape[0])
                self.laminar_flame_velocity = self.laminar_flame_velocity.append(add_row)

        # Calculation with Cantera
        else:

            # Define the fuel from the surrogate composition in a ARCANE friendly style
            fuel = "X/"
            for spec in self.composition:
                fuel += f"{spec}/{format(self.composition[spec], '.4f')}/"
            fuel = fuel[:-1]

            # Define cases to run
            caselist = []
            caselist.extend(cases.create_case(reactor="C1DP",
                                              mechanism=self.mechanism,
                                              fuel=fuel,
                                              oxidizer=oxidizer,
                                              pressure=str(pressure),
                                              temperature=str(temperature),
                                              phi=str(phi)))

            cases.run_cases(caselist, self.mechanism, overwrite=False)

            # Store the data
            for case in caselist:
                data_dict_from_case = case.data_dict(self.mechanism)
                add_row = pd.Series(data=[data_dict_from_case['Temperature'][0],
                                          data_dict_from_case['Pressure'][0],
                                          case.parameters['phi'],
                                          case.extract_quantity(self.mechanism, 'Laminar flame speed')],
                                    index=self.laminar_flame_velocity.columns,
                                    name=self.laminar_flame_velocity.shape[0])
                self.laminar_flame_velocity = self.laminar_flame_velocity.append(add_row)

        return

    def calc_ignition_delay_time_from_target(self, target, oxidizer="X/O2/1.0/N2/14.28"):
        """
        Calculate the ignition delay time with Cantera based on the thermodynamic conditions defined for the target

        Parameters
        ----------

        target:
            Surrogate object to match
        oxidizer:
            oxidizer string in ARCANE case definition style

        Created: 21/07/16 [TL]
        """
        # Erase data in the dataframe
        self.ignition_delay_time.drop(self.ignition_delay_time.index, inplace=True)

        # Define the fuel from the surrogate composition and oxidizer
        fuel = "X/"
        for spec in self.composition:
            fuel += f"{spec}/{format(self.composition[spec], '.4f')}/"
        fuel = fuel[:-1]

        # Formatting thermo conditions to ARCANE cases definition
        added_pressure = []
        added_temperature = []
        added_phi = []
        pressure = ""
        temperature = ""
        phi = ""
        for i in range(len(target.ignition_delay_time['temperature'])):
            if not target.ignition_delay_time['temperature'][i] in added_temperature:
                added_temperature.append(target.ignition_delay_time['temperature'][i])
                temperature += f'{target.ignition_delay_time["temperature"][i]}-'
            if not target.ignition_delay_time['pressure'][i] in added_pressure:
                added_pressure.append(target.ignition_delay_time['pressure'][i])
                pressure += f'{target.ignition_delay_time["pressure"][i]}-'
            if not target.ignition_delay_time['phi'][i] in added_phi:
                added_phi.append(target.ignition_delay_time['phi'][i])
                phi += f'{target.ignition_delay_time["phi"][i]}-'
        phi = phi[:-1]
        temperature = temperature[:-1]
        pressure = pressure[:-1]

        # Define cases to run
        caselist = []
        caselist.extend(cases.create_case(reactor="C0DV",
                                          mechanism=self.mechanism,
                                          fuel=fuel,
                                          oxidizer=oxidizer,
                                          pressure=str(pressure),
                                          temperature=str(temperature),
                                          phi=str(phi)))

        cases.run_cases(caselist, self.mechanism, overwrite=False)

        # Store the data
        for case in caselist:
            data_dict_from_case = case.data_dict(self.mechanism)
            add_row = pd.Series(data=[data_dict_from_case['Temperature'][0],
                                      data_dict_from_case['Pressure'][0],
                                      case.parameters['phi'],
                                      case.extract_quantity(self.mechanism, 'Ignition delay time')],
                                index=self.ignition_delay_time.columns,
                                name=self.ignition_delay_time.shape[0])
            self.ignition_delay_time = self.ignition_delay_time.append(add_row)

        return

    def calc_ignition_delay_time(self, temperature, pressure, phi, oxidizer="X/O2/1.0/N2/14.28"):
        """
        Calculate the ignition delay time with Cantera with the thermodynamic conditions in arguments

        Parameters
        ----------

        temperature:
            temperature string in ARCANE case definition style
        pressure:
            pressure string in ARCANE case definition style
        phi:
            phi string in ARCANE case definition style
        oxidizer:
            oxidizer string in ARCANE case definition style

        Created: 21/07/16 [TL]
        """

        # Define the fuel from the surrogate composition and oxidizer
        fuel = "X/"
        for spec in self.composition:
            fuel += f"{spec}/{format(self.composition[spec], '.4f')}/"
        fuel = fuel[:-1]

        # Define cases to run
        caselist = []
        caselist.extend(cases.create_case(reactor="C0DV",
                                          mechanism=self.mechanism,
                                          fuel=fuel,
                                          oxidizer=oxidizer,
                                          pressure=str(pressure),
                                          temperature=str(temperature),
                                          phi=str(phi)))

        cases.run_cases(caselist, self.mechanism, overwrite=False)

        # Store the data
        for case in caselist:
            data_dict_from_case = case.data_dict(self.mechanism)
            add_row = pd.Series(data=[data_dict_from_case['Temperature'][0],
                                      data_dict_from_case['Pressure'][0],
                                      case.parameters['phi'],
                                      case.extract_quantity(self.mechanism, 'Ignition delay time')],
                                index=self.ignition_delay_time.columns,
                                name=self.ignition_delay_time.shape[0])
            self.ignition_delay_time = self.ignition_delay_time.append(add_row)

        return

    def add_ignition_delay_time(self, temperature, pressure, phi, tig):
        """
        Adds external data of ignition delay time for a surrogate

        Parameters
        ----------

        temperature:
            list of temperatures
        pressure:
            list of pressure
        phi:
            list of stoichiometric ratio
        tig:
            list of the corresponding ignition delay time

        Created: 21/07/16 [TL]
        """

        # Check if the data are well formatted
        if isinstance(temperature, list):
            if len(temperature) != len(pressure) or len(temperature) != len(phi) or len(temperature) != len(tig):
                logger.error("The arguments of add_ignition_delay_time should be lists of the same length")
                quit()

            # Store the data
            for elem, temp in enumerate(temperature):
                add_row = pd.Series(data=[temp,
                                          pressure[elem],
                                          phi[elem],
                                          tig[elem]],
                                    index=self.ignition_delay_time.columns,
                                    name=self.ignition_delay_time.shape[0])
                self.ignition_delay_time = self.ignition_delay_time.append(add_row)
        # else:
        #     add_row = pd.Series(data=[temperature,
        #                               pressure,
        #                               phi,
        #                               tig],
        #                         index=self.ignition_delay_time.columns,
        #                         name=self.ignition_delay_time.shape[0])
        #     self.ignition_delay_time = self.ignition_delay_time.append(add_row)

        return

    def add_laminar_flame_velocity(self, temperature, pressure, phi, sl):
        """
        Adds external data of laminar flame speed for a surrogate

        Parameters
        ----------

        temperature:
            list of temperatures
        pressure:
            list of pressure
        phi:
            list of stoichiometric ratio
        sl:
            list of the corresponding laminar flame speed

        Created: 21/07/16 [TL]
        """
        # Check if the data are well formatted
        if isinstance(temperature, list):
            if len(temperature) != len(pressure) or len(temperature) != len(phi) or len(temperature) != len(sl):
                logger.error("The arguments of add_laminar_flame_velocity should be lists of the same length")
                quit()

            # Store the data
            for elem, temp in enumerate(temperature):
                add_row = pd.Series(data=[temp,
                                          pressure[elem],
                                          phi[elem],
                                          sl[elem]],
                                    index=self.laminar_flame_velocity.columns,
                                    name=self.laminar_flame_velocity.shape[0])
                self.laminar_flame_velocity = self.laminar_flame_velocity.append(add_row)
        # else:
        #     add_row = pd.Series(data=[temperature,
        #                               pressure,
        #                               phi,
        #                               sl],
        #                         index=self.laminar_flame_velocity.columns,
        #                         name=self.laminar_flame_velocity.shape[0])
        #     self.laminar_flame_velocity = self.laminar_flame_velocity.append(add_row)

        return


    def calc_distillation_curve(self, init_guess=450, mixture_method="ideal"):
        """
        Add a distillation curve as a property of the Surrogate class,
        returns list of two lists in the key 'distillation_curve' of the Surrogate.properties dictionary,
        first list : vapor fraction evaporated, second list : boiling temperature corresponding

        Parameters
        ----------

        init_guess:
            temperature guess for the boiling temperature at 1% volume evaporated

        Created: 21/07/01 [TL]
        """
        # Test if the composition is given by species
        if self.raw_composition:
            logger.error('Composition should given by species to perform distillation curve calculation')
            quit()

        # Test if the vapor pressure is already given or if it is known by ARCANE
        for spec in self.composition:
            if not 'vapor_pressure' in self.species_properties[spec]:
                if self.check_chem_db:
                    self.species_properties[spec]['vapor_pressure'] = AVBP.species_sat_pressure(
                        self.species_properties[spec])
                else:
                    logger.error('vapor_pressure should be given for the species')
                    quit()
                       
        # Test of method
        if mixture_method != "ideal" and mixture_method != "non-ideal":
            logger.error("Mixture method for distillation curve must be 'ideal' or 'non-ideal'")
            quit()

        # Local variable initialization
        number_species = len(self.composition.keys())
        x_l = np.zeros([number_species, 100])
        x_v = np.zeros([number_species, 99])
        boiling_temperature = np.zeros(99)
        vap_frac = np.zeros(99)
        fuel = list(self.composition.keys())
        volume = 1
        ratio_p = np.zeros(number_species)
        ptot = 101325  # standard distillation curve are given for 1 atm
        tpf_saturation_pressure = []

        # Initialisation for activity coefficient calculation if required
        chemgroups = []
        if mixture_method == "non-ideal":
            # Test if thermo module is available
            try :
                from thermo.unifac import UFIP, UFSG, UNIFAC
            except ModuleNotFoundError :
                logger.error("Module thermo must be installed in your python environment for non-ideal mixture ! (pip install thermo)")
                quit()
            for name in fuel:
                chemgroups.append(species_chemical_groups(name))

        # Formatting each species vapor pressure into a concatenated list
        # Liquid molar fraction initialization
        i = 0
        for name in fuel:
            tpf_saturation_pressure.append(self.species_properties[name]['vapor_pressure'])
            x_l[i, 0] = self.composition[name]
            i = i + 1

        # Loop for iterative calculation of the boiling_temperature resulting of the liquid-vapor equilibrium
        i = 0
        while i < 99:

            # Get boiling_temperature from the Raoult's equation
            if i == 0:
                if mixture_method == 'non-ideal':
                    boiling_temperature[i] = fsolve(non_ideal_raoult_law, init_guess, args=(tpf_saturation_pressure, x_l[:, 0], ptot, chemgroups))
                else :
                    boiling_temperature[i] = fsolve(raoult_law, init_guess, args=(tpf_saturation_pressure, x_l[:, 0], ptot))
            else:
                if mixture_method == 'non-ideal':
                    boiling_temperature[i] = fsolve(non_ideal_raoult_law, boiling_temperature[i - 1],
                                                args=(tpf_saturation_pressure, x_l[:, i], ptot, chemgroups))
                else :
                    boiling_temperature[i] = fsolve(raoult_law, boiling_temperature[i - 1],
                                                args=(tpf_saturation_pressure, x_l[:, i], ptot))

            # Defined the evaporated fraction associated with boiling_temperature
            vap_frac[i] = i + 0.5

            # Solve for each vapor molar fraction to update the liquid molar fraction
            for j in range(len(tpf_saturation_pressure)):
                ratio_p[j] = psat(tpf_saturation_pressure[j], boiling_temperature[i]) / ptot
            x_v[:, i] = x_l[:, i] * ratio_p[:]

            # Update the liquid molar fraction considering that at the step 1% of the initial volume as evaporated
            x_l[:, i + 1] = (volume * x_l[:, i] - 0.01 * x_v[:, i]) / (volume - 0.01)
            x_l[:, i + 1] = x_l[:, i + 1] / sum(x_l[:, i + 1])

            # Update loop variables
            volume = volume - 0.01
            i = i + 1

        # Set the distillation curve as a surrogate property
        self.properties['distillation_curve'] = [vap_frac, boiling_temperature - 273.15]

        return


def psat(tpf_saturation_pressure, temperature):
    """
    Function calculating the vapor pressure for a given temperature

    Parameters
    ----------

    tpf_saturation_pressure:
        vapor pressure list (from 0K to critical temperature by step of 10K, AVBP style)
    temperature:
        temperature

    """
    ind = int(0.1 * temperature) - 1 # Apparently there is a shift compared to AVBP (validated with http://www.vle-calc.com/distillation_calculation.html)
    fr = temperature * 0.1 - float(ind) - 1.0
    try:
        vapor_pressure = tpf_saturation_pressure[ind + 1] * (1.0 - fr) + tpf_saturation_pressure[ind + 2] * fr
    except:
        logger.warning("One species has got over its critical pressure (psat is set to pcrit)")
        vapor_pressure = tpf_saturation_pressure[-2]
    return vapor_pressure


def raoult_law(temperature, tpf_saturation_pressure, x_l, gas_pressure):
    """
    Function calculating the Raoult's law, at equilibrium returns zero

    Parameters
    ----------

    temperature:
        temperature
    tpf_saturation_pressure:
        list of vapor pressure lists for each species
    x_l:
        list of species molar fraction for each species
    gas_pressure:
        ambient pressure

    """
    f = -1
    for i in range(len(tpf_saturation_pressure)):
        f = f + x_l[i] / gas_pressure * psat(tpf_saturation_pressure[i], temperature)
    return f


def non_ideal_raoult_law(temperature, tpf_saturation_pressure, x_l, gas_pressure, chemgroups):
    """
    Function calculating the Raoult's law, at equilibrium returns zero

    Parameters
    ----------

    temperature:
        temperature
    tpf_saturation_pressure:
        list of vapor pressure lists for each species
    x_l:
        list of species molar fraction for each species
    gas_pressure:
        ambient pressure
    chemgroups:
        chemical groups for UNIFAC activity coefficient calculation 

    """
    from thermo.unifac import UFIP, UFSG, UNIFAC
    GE = UNIFAC.from_subgroups(chemgroups=chemgroups, T = temperature, xs=x_l, version=0, interaction_data=UFIP, subgroups=UFSG)
    activity_coef = GE.gammas()

    f = -1
    for i in range(len(tpf_saturation_pressure)):
        f = f + activity_coef[i] * x_l[i] / gas_pressure * psat(tpf_saturation_pressure[i], temperature)
    return f


    # def get_potential_palette_species(self):
    #     """
    #     Identifies the hydrocarbon species in the mechanism that belong to a specific family
    #
    #
    #     :return: list of potential species that can be used for the surrogate generation
    #     """
    #
    #     for spec in self.mechanism.species_names:
    #         carbon_atoms =
    #         hydrogen_atoms =

    # def hydrogen_carbon_ratio(self):

def species_chemical_groups(species_name):
    """
    Temporary store of chemical groups for liquid species for activity coeffcient calculation
    """
    # 1 : CH3
    # 2 : CH2
    # 3 : CH
    # 4 : C
    # 9 : ACH (Aromatics)
    # 10 : AC
    chemical_groups_dict={
        'NC12H26' : {1:2, 2:10},
        'MCYC6' : {1:1, 2:5, 3:1},
        'TETRALIN' : {2:4, 9:4, 10:2},
        "C7H8": {1:1, 9:5, 10:1},
        "NC7H16" : {1:2, 2:5},
        "XYLENE" : {1:2, 9:4, 10:2},
        "IC8H18" : {1:5, 2:1, 3:1, 4:1},
        "DECALIN" : {2:8, 3:2},
        "NC10H22" : {1:2, 2:8},
        "C11H10" : {1:1, 9:7, 10:3},
        "IC12H26" : {1:7, 2:2, 3:1, 4:2},
        "IC16H34" : {1:9, 2:3, 3:1,4:3},
        "NC16H34" : {1:2,2:14},
        "TMBENZ" : {1:4, 9:2, 10:4},
        "C10H7CH3" : {1:1, 9:7, 10:3}
    }

    if not species_name in chemical_groups_dict.keys():
        logger.error("The species does not exist in the chemical_groups_dict of surrogate.py, Please add it")
        quit()

    return chemical_groups_dict[species_name]