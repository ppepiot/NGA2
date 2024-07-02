import ARCANE.display as display

logger = display.Logger()
logger.set_log('logCompute')


# (v)('m')(v)
# Look for the Zoidberg for customisable parameters
# Feel free to search it and customise your output


def write(ctmech, output_name=None, transport=True, thermo=True):
    """Writes the .inp file containing the kinetics in Chemkin format and associated transport and thermodynamic
    files if specified

    Parameters
    ----------
    ctmech :
        class `Cantera.Solution` object
    output_name : str
        name of the output mechanism (Default value = None)
    transport : bool
        if True, writes a .tran file containing the transport data (Default value = True)
    thermo : bool
        if True, writes a .thermo file containing the thermodynamic data (Default value = True)

    """

    if not output_name:
        output_name = "chemkin.inp"

    file = open(output_name, 'w+')

    # Writing the elements
    elements_list = ctmech.element_names

    elements_string = "ELEMENTS\n"
    for element in elements_list:
        elements_string += f"{element}\n"
    elements_string += "END\n\n"

    file.write(elements_string.upper())

    # Writing the species
    species_list = ctmech.species_names
    # Format parameters
    number_of_columns = 3  # (v)('m')(v)
    spec_length = 30  # (v)('m')(v)

    counter = 0
    species_string = "SPECIES\n"
    for spec in species_list:
        species_string += f"{spec:{spec_length}}"
        counter += 1
        if (counter % number_of_columns) == 0:
            species_string += '\n'
        elif spec == species_list[-1]:
            species_string += '\n'
    species_string += '\nEND\n\n'

    file.write(species_string)

    # Write reactions
    reactions_string = "REACTIONS\n\n"
    file.write(reactions_string)

    # Retrieving reactions
    reactions = ctmech.reactions()

    # Length parameters
    equation_length = 52  # (v)('m')(v)
    pre_exp_length = 11  # (v)('m')(v)
    beta_length = 11  # (v)('m')(v)
    energy_length = 11  # (v)('m')(v)
    parameter_format = ".3E"  # (v)('m')(v)
    efficiency_format = ".2f"  # (v)('m')(v)
    orders_format = ".2f"  # (v)('m')(v)

    calories_constant = 4184.0  # number of calories in 1000 Joules of energy

    for reaction in reactions:
        reaction_type = reaction.reaction_type

        reaction_equation = reaction.equation
        reaction_equation = reaction_equation.replace(' ', '').replace('<=>', '=')

        # Determining the reaction order to apply the correct rate correction
        if reaction.orders:
            order = sum(reaction.orders.values())
        else:
            order = sum(reaction.reactants.values())

        if reaction_type == 2:
            order += 1

        # Initialised for duplicity
        duplicate_string = ""

        # Initialised for arbitrary orders
        orders_string = ""

        # Initialised for type 2 and 4
        efficiency_string = ""

        # Initialised for type 4
        low_rate_string = ""
        falloff_string = ""

        # Initialised for type 5
        plog_string = ""

        # Elementary reaction
        if reaction_type in [1, 2]:

            rate = reaction.rate
            pre_exp = rate.pre_exponential_factor * 10 ** ((order - 1) * 3)
            beta = rate.temperature_exponent
            energy = rate.activation_energy / calories_constant

        elif reaction_type == 4:

            # Investigate the correct physical meaning
            products_order = sum(reaction.products.values())

            # High rates
            high_rate = reaction.high_rate

            high_A = high_rate.pre_exponential_factor * 10 ** ((order - 1) * 3)
            high_b = high_rate.temperature_exponent
            high_Ea = high_rate.activation_energy / calories_constant

            # Low rates
            low_rate = reaction.low_rate

            low_A = low_rate.pre_exponential_factor * 10 ** (order * 3)
            low_b = low_rate.temperature_exponent
            low_Ea = low_rate.activation_energy / calories_constant

            pre_exp = high_A
            beta = high_b
            energy = high_Ea

            low_rate_string = f" LOW/ {low_A:{pre_exp_length}{parameter_format}}" \
                              f"{low_b:{beta_length}{parameter_format}}{low_Ea:>{energy_length}{parameter_format}}/"

            falloff_type = reaction.falloff.type

            falloff_length = 11  # (v)('m')(v)
            falloff_format = ".4E"  # (v)('m')(v)

            if falloff_type == 'Troe':

                falloff_coefficients = reaction.falloff.parameters

                if falloff_coefficients[3] == 0:
                    falloff_string = f"TROE/ {falloff_coefficients[0]:{falloff_length}{falloff_format}}" \
                                     f"{falloff_coefficients[1]:{falloff_length}{falloff_format}}" \
                                     f"{falloff_coefficients[2]:{falloff_length}{falloff_format}}/"
                else:
                    falloff_string = f"TROE/ {falloff_coefficients[0]:{falloff_length}{falloff_format}}" \
                                     f"{falloff_coefficients[1]:{falloff_length}{falloff_format}}" \
                                     f"{falloff_coefficients[2]:{falloff_length}{falloff_format}}" \
                                     f"{falloff_coefficients[3]:{falloff_length}{falloff_format}}/"

            elif falloff_type == 'Simple':
                falloff_string = ""

            elif falloff_type == 'SRI':

                falloff_coefficients = reaction.falloff.parameters

                falloff_string = f"SRI/ {falloff_coefficients[0]:{falloff_length}{falloff_format}}" \
                                 f"{falloff_coefficients[1]:{falloff_length}{falloff_format}}" \
                                 f"{falloff_coefficients[2]:{falloff_length}{falloff_format}}"\
                                 f"{falloff_coefficients[3]:{falloff_length}{falloff_format}}" \
                                 f"{falloff_coefficients[4]:{falloff_length}{falloff_format}}/"

            else:
                logger.error(f'This type of falloff parametrisation ({falloff_type}) is not coded.')

        elif reaction_type == 5:
            plog_rates = reaction.rates

            pressure_length = 20  # (v)('m')(v)
            pressure_format = ".6E"  # (v)('m')(v)

            for rate in plog_rates:
                pressure = rate[0] / 1.01325e5  # Pascal to atm
                pre_exp_plog = rate[1].pre_exponential_factor * 10 ** ((order - 1) * 3)
                beta_plog = rate[1].temperature_exponent
                energy_plog = rate[1].activation_energy / calories_constant

                plog_string += f" PLOG/ {pressure:{pressure_length}{pressure_format}}" \
                               f"{pre_exp_plog:{pre_exp_length}{parameter_format}}" \
                               f"{beta_plog:{beta_length}{parameter_format}}" \
                               f"{energy_plog:>{energy_length}{parameter_format}}/\n"

                # Rate coefficients following the equation taken at 1 atm
                if pressure == 1:
                    pre_exp = pre_exp_plog
                    beta = beta_plog
                    energy = energy_plog

        else:
            logger.error("The converter encountered an unknown reaction type "
                         "that needs to be coded or manually written.")

        # Efficiencies for type 2 and 4
        if reaction_type in [2, 4]:
            efficiencies = reaction.efficiencies
            for spec in efficiencies:
                efficiency_string += f"{spec}/ {efficiencies[spec]:{efficiency_format}} "

        # Duplicity of the reaction
        if reaction.duplicate:
            duplicate_string = " DUPLICATE"

        # Arbitrary orders
        if reaction.orders:
            orders_string = " FORD/"
            for spec in reaction.orders:
                orders_string += f"{spec} {reaction.orders[spec]:{orders_format}}/"

        reactions_string = f"{reaction_equation:{equation_length}}"
        reactions_string += f"{pre_exp:{pre_exp_length}{parameter_format}}{beta:{beta_length}{parameter_format}}" \
                            f"{energy:>{energy_length}{parameter_format}}"

        # Writing equation and rate coefficents
        file.write(reactions_string)
        file.write("\n")
        # Writing duplicity
        if duplicate_string:
            file.write(duplicate_string)
            file.write("\n")
        # Writing arbitrary orders
        if orders_string:
            file.write(orders_string)
            file.write("\n")
        # Writing low rate coefficients for type 4
        if low_rate_string:
            file.write(low_rate_string)
            file.write("\n")
        # Writing falloff parameters for type 4
        if falloff_string:
            file.write(falloff_string)
            file.write("\n")
        # Writing efficiency if type 2 or 4
        if efficiency_string:
            file.write(efficiency_string)
            file.write("\n")
        # Writing plog rates if type 5
        if plog_string:
            file.write(plog_string)

        file.write("\n")

    # Writing end statement
    file.write("END")

    logger.goodnews(f"Kinetics file correctly written as {output_name}")

    if thermo:
        thermo_name = output_name.replace('.inp', '.thermo')
        write_thermo(ctmech, output_name=thermo_name)

    if transport:
        transport_name = output_name.replace('.inp', '.tran')
        write_transport(ctmech, output_name=transport_name)


def write_transport(ctmech, output_name=None):
    """Writes the .tran file containing the transport data in Chemkin format

    Parameters
    ----------
    ctmech :
        class `Cantera.Solution` object
    output_name : str
        name of the output mechanism (Default value = None)

    """
    if ctmech.transport_model == 'Transport':
        logger.info("No transport file will be written has the cti does not contain transport data.")

    if not output_name:
        output_name = "chemkin.tran"

    file = open(output_name, 'w+')

    species = ctmech.species()

    for spec in species:
        name = spec.name

        geometry = spec.transport.geometry
        # Checking geometry
        if geometry == 'atom':
            geometry_number = 0
        elif geometry == 'linear':
            geometry_number = 1
        elif geometry == 'nonlinear':
            geometry_number = 2

        # Retrieving the parameters and scaling them
        well_depth = spec.transport.well_depth / 1.380649e-23  # Joules to Kelvin (Bolzmann constant = 1.380649e-23)
        diameter = spec.transport.diameter * 1e10  # meters to Angström
        dipole_moment = spec.transport.dipole / (1e-21 / 299792458)  # Coulom.meters to Debyes (c = 299792458 m/s)
        polarisation = spec.transport.polarizability * 1e30  # cubic meters to cubic Angström
        rotational_relax = spec.transport.rotational_relaxation

        # Customisable parameters
        name_length = 15  # (v)('m')(v)
        geometry_length = 10  # (v)('m')(v)
        values_format = '.3f'  # (v)('m')(v)
        well_depth_length = 10  # (v)('m')(v)
        diameter_length = 10  # (v)('m')(v)
        dipole_moment_length = 10  # (v)('m')(v)
        polarisation_length = 10  # (v)('m')(v)
        rotational_relax_length = 10  # (v)('m')(v)

        transport_string = f"{name:{name_length}}{geometry_number:{geometry_length}}" \
                           f"{well_depth:{well_depth_length}{values_format}}" \
                           f"{diameter:{diameter_length}{values_format}}" \
                           f"{dipole_moment:{dipole_moment_length}{values_format}}" \
                           f"{polarisation:{polarisation_length}{values_format}}" \
                           f"{rotational_relax:{rotational_relax_length}{values_format}}\n"

        file.write(transport_string)

    logger.goodnews(f"Transport file correctly written as {output_name}")


def write_thermo(ctmech, output_name=None):
    """Writes the .thermo file containing the thermodynamic data in Chemkin format

    Parameters
    ----------
    ctmech :
        class `Cantera.Solution` object
    output_name : str
        name of the output mechanism (Default value = None)

    """

    if not output_name:
        output_name = "chemkin.thermo"

    file = open(output_name, 'w+')

    # Writing header
    thermo_string = "THERMO ALL\n"

    # Customisable parameters
    min_temperature = 270.  # (v)('m')(v)
    intermediate_temperature = 1000.  # (v)('m')(v)
    max_temperature = 5000.  # (v)('m')(v)
    temperatures_length = 10
    temperatures_format = '.1f'

    thermo_string += f"{min_temperature:{temperatures_length}{temperatures_format}}" \
                     f"{intermediate_temperature:{temperatures_length}{temperatures_format}}" \
                     f"{max_temperature:{temperatures_length}{temperatures_format}}\n"

    file.write(thermo_string)

    species = ctmech.species()

    # Customisable parameters
    species_length = 20

    for spec in species:
        species_string = f"{spec.name:{species_length}}"

        # Creating the composition string
        composition = spec.composition
        composition_string = ""
        element_length = 2  # (v)('m')(v)
        element_number_length = 3  # (v)('m')(v)
        for element in composition:
            element_number = int(composition[element])
            composition_string += f"{element:{element_length}}{element_number:>{element_number_length}}"

        species_string += composition_string.upper()

        # Adding species phase
        phase = " G"
        phase_length = 30 - len(composition_string)  # (v)('m')(v)

        species_string += f"{phase:>{phase_length}}"

        # Temperatures boundaries
        nasa_coefficients = spec.thermo.coeffs

        min_temperature = spec.thermo.min_temp
        max_temperature = spec.thermo.max_temp
        intermediate_temperature = nasa_coefficients[0]

        species_string += f"{min_temperature:{temperatures_length}{temperatures_format}}" \
                          f"{intermediate_temperature:{temperatures_length}{temperatures_format}}" \
                          f"{max_temperature:{temperatures_length}{temperatures_format}} 1\n"

        # NASA polynomials
        nasa_coeffs_length = 16
        nasa_coeffs_format = '.8E'
        second_row_string = f"{nasa_coefficients[1]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                            f"{nasa_coefficients[2]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                            f"{nasa_coefficients[3]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                            f"{nasa_coefficients[4]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                            f"{nasa_coefficients[5]:{nasa_coeffs_length}{nasa_coeffs_format}} 2\n"

        third_row_string = f"{nasa_coefficients[6]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                           f"{nasa_coefficients[7]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                           f"{nasa_coefficients[8]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                           f"{nasa_coefficients[9]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                           f"{nasa_coefficients[10]:{nasa_coeffs_length}{nasa_coeffs_format}} 3\n"

        if len(nasa_coefficients) == 15:
            filling = " " * nasa_coeffs_length
            forth_row_string = f"{nasa_coefficients[11]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[12]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[13]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[14]:{nasa_coeffs_length}{nasa_coeffs_format}}{filling} 4\n"

            fifth_row_string = ""
        else:
            forth_row_string = f"{nasa_coefficients[11]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[12]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[13]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[14]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[15]:{nasa_coeffs_length}{nasa_coeffs_format}} 4\n"

            filling = " " * 2 * nasa_coeffs_length
            fifth_row_string = f"{nasa_coefficients[16]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[17]:{nasa_coeffs_length}{nasa_coeffs_format}}" \
                               f"{nasa_coefficients[18]:{nasa_coeffs_length}{nasa_coeffs_format}}{filling} 5\n"

            logger.warning("This script was not designed for NASA9 polynomials, some adjustments must be made...")

        species_string += second_row_string + third_row_string + forth_row_string + fifth_row_string
        file.write(species_string)

    # Ending
    file.write("END")

    logger.goodnews(f"Thermodynamic file correctly written as {output_name}")
