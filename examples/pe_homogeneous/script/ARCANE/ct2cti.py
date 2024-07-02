"""Script to write cti files from Parker Clayton
Source: https://groups.google.com/forum/#!topic/cantera-users/67860KYm6JY/discussion

TODO
----
- Figure out how to best use/redo/give credit to this script
"""

from __future__ import division, print_function

import textwrap
from string import Template

import cantera as ct


def write(solution, where='cantera.cti', species_change_dict={}, species_qss=[], dummy='', kinetics=None, transport=None):
    """Function to write cantera solution object to cti file.

    Parameters
    ----------
    solution :
        class `Cantera.Solution` object
    where : str
        path to the written file (Default value = 'cantera.cti')
    species_change_dict : dict
        dictionary of species names that needs to be replaced (Default value = {})
    species_qss : list of str
        QSS species that will be written in a comment (Default value = [])
    dummy : str
        fake reaction for Cantera (Default value = '')
    kinetics : str
        kinetic model (Default value = None)
    transport : str
        transport model

    Returns
    -------
    output_file_name : str
        Name of trimmed mechanism file (.cti) (Default value = None)

    """

    trimmed_solution = solution
    input_file_name_stripped = trimmed_solution.name
    if not input_file_name_stripped:
        input_file_name_stripped = "gas"

    if not transport:
        transport_model = solution.transport_model
    else:
        transport_model = transport

    if not where.endswith('.cti'):
        where += '.cti'
    output_file_name = where

    with open(output_file_name, 'w+') as f:

        # Get solution temperature and pressure
        solution_temperature = '{:.6E}'.format(300.0, 1)
        solution_pressure = '{:.6E}'.format(ct.one_atm, 1)

        # Work Functions

        # number of calories in 1000 Joules of energy
        calories_constant = 4184.0

        def eliminate(input_string, char_to_replace, spaces='single'):
            """Eliminate characters from a string

            Parameters
            ----------
            input_string :
                string to be modified
            char_to_replace :
                array of character strings to be removed
            spaces :
                string to know if one or double spaced (Default value = 'single')

            Returns
            -------
            output_string :
                modified string

            """
            for char in char_to_replace:
                output_string = input_string.replace(char, "")
            if spaces == 'double':
                output_string = input_string.replace(" ", "  ")
            return output_string

        def wrap_nasa(input_string):
            """Wrap string to cti NASA format width

            Parameters
            ----------
            input_string :
                string to be modified

            Returns
            -------
            output_string:
                string modified

            """
            output_string = textwrap.fill(
                input_string,
                width=50,
                subsequent_indent=16 * ' '
            )
            return output_string

        def section_break(title):
            """Insert break and new section title into cti file

            Parameters
            ----------
            title :
                title string for next section_break

            """
            f.write('#' + '-' * 75 + '\n')
            f.write('#  ' + title + '\n')
            f.write('#' + '-' * 75 + '\n\n')

        def replace_multiple(input_string, replace_list):
            """Replace multiple characters in a string

            Parameters
            ----------
            input_string :
                string to be modified
            replace_list :
                list containing items to be replaced (value replaces key)

            Returns
            -------
            input_string :
                string modified

            """
            for original_character, new_character in replace_list.items():
                input_string = input_string.replace(original_character,
                                                    new_character)
            return input_string

        def build_arrhenius(equation_object, equation_type):
            """Builds Arrhenius coefficient string

            Parameters
            ----------
            equation_object :
                `Cantera` equation object
            equation_type :
                string of equation type, either 'ElementaryReaction', 'ThreeBodyReaction'

            Returns
            -------
            Arrhenius coefficient as string

            """
            coeff_sum = sum(equation_object.reactants.values())
            pre_exponential_factor = equation_object.rate.pre_exponential_factor
            temperature_exponent = equation_object.rate.temperature_exponent
            activation_energy = equation_object.rate.activation_energy / calories_constant

            if equation_type == 'ElementaryReaction':
                if coeff_sum == 1:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor))
                if coeff_sum == 2:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 3))
                if coeff_sum == 3:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 6))
            if equation_type == 'ThreeBodyReaction':
                if coeff_sum == 1:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 3))
                if coeff_sum == 2:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 6))

            if (equation_type != 'ElementaryReaction'
                    and equation_type != 'ThreeBodyReaction'):
                pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor))

            temperature_exponent = str('{:.6E}'.format(temperature_exponent))
            activation_energy = str('{:.6E}'.format(activation_energy))

            arrhenius = [pre_exponential_factor,
                         temperature_exponent,
                         activation_energy]

            return str(arrhenius).replace("\'", "")

        def build_modified_arrhenius(equation_object, t_range):
            """Builds Arrhenius coefficient strings for high and low temperature ranges

            Parameters
            ----------
            equation_object :
                `Cantera` equation object
            t_range :
                simple string ('high' or 'low') to designate temperature range

            Returns
            -------
            Modified Arrhenius coefficient as string

            """

            if t_range == 'high':
                pre_exponential_factor = equation_object.high_rate.pre_exponential_factor
                temperature_exponent = equation_object.high_rate.temperature_exponent
                activation_energy = equation_object.high_rate.activation_energy / calories_constant

                if sum(equation_object.reactants.values()) == 1:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor))
                elif sum(equation_object.reactants.values()) == 2:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 3))
                elif sum(equation_object.reactants.values()) == 3:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 6))
                else:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor))

                temperature_exponent = str('{:.6E}'.format(temperature_exponent))
                activation_energy = str('{:.6E}'.format(activation_energy))

                arrhenius_high = [pre_exponential_factor,
                                  temperature_exponent,
                                  activation_energy]

                return str(arrhenius_high).replace("\'", "")

            elif t_range == 'low':
                pre_exponential_factor = equation_object.low_rate.pre_exponential_factor
                temperature_exponent = equation_object.low_rate.temperature_exponent
                activation_energy = equation_object.low_rate.activation_energy / calories_constant

                if sum(equation_object.reactants.values()) == 1:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 3))
                elif sum(equation_object.reactants.values()) == 2:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 6))
                elif sum(equation_object.reactants.values()) == 3:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor * 10 ** 9))
                else:
                    pre_exponential_factor = str('{:.6E}'.format(pre_exponential_factor))

                temperature_exponent = str('{:.6E}'.format(temperature_exponent))
                activation_energy = str('{:.6E}'.format(activation_energy))

                arrhenius_low = [pre_exponential_factor,
                                 temperature_exponent,
                                 activation_energy]

                return str(arrhenius_low).replace("\'", "")

            else:
                return None

        def build_falloff(equation_object):
            """Creates falloff reaction Troe parameter string

            Parameters
            ----------
            equation_object :
                `Cantera` reaction object

            Returns
            -------
            falloff_string :
                string associated to fall-out

            """

            if equation_object.falloff.type == 'Simple':
                falloff_string = ')\n\n'
            elif equation_object.falloff.type == "Troe":
                parameters = equation_object.falloff.parameters
                falloff_string = str(
                    ',\n        falloff = Troe(' +
                    'A = ' + str(parameters[0]) +
                    ', T3 = ' + str(parameters[1]) +
                    ', T1 = ' + str(parameters[2]) +
                    ', T2 = ' + str(parameters[3]) + ')       )\n\n'
                )
            elif equation_object.falloff.type == "SRI":
                parameters = equation_object.falloff.parameters
                falloff_string = str(
                    ',\n        falloff = SRI(' +
                    'A = ' + str(parameters[0]) +
                    ', B = ' + str(parameters[1]) +
                    ', C = ' + str(parameters[2]) +
                    ', D = ' + str(parameters[3]) +
                    ', E = ' + str(parameters[4]) + ')       )\n\n'
                )
            else:
                print("WARNING: This Falloff type is not coded yet")

            if species_qss:
                falloff_string = falloff_string.replace('        falloff', '#        falloff')

            return falloff_string

        def build_species_string(species_names, commented=False):
            """Formats species list at top of mechanism file

            Parameters
            ----------
            species_names :
                string names of the species
            commented :
                write comments or not (Default value = False)

            Returns
            -------
            output_file_name:
                name of the output file

            """
            species_list_string = ''
            line = 1

            if commented:
                suffix = '#'
            else:
                suffix = ''

            for sp_str in species_names:
                # get length of string next species is added
                length_new = len(sp_str)
                length_string = len(species_list_string)
                total = length_new + length_string + 3
                # if string will go over width, wrap to new line
                if line == 1:
                    if total >= 55:
                        species_list_string += '\n'
                        species_list_string += suffix + ' ' * 17
                        line += 1
                if line > 1:
                    if total >= 70 * line:
                        species_list_string += '\n'
                        species_list_string += suffix + ' ' * 17
                        line += 1
                species_list_string += sp_str + ' '
            return species_list_string

        # Write title block to file

        section_break('CTI File converted from Solution Object')
        unit_string = "units(length = \"cm\", time = \"s\"," + \
                      " quantity = \"mol\", act_energy = \"cal/mol\")"
        f.write(unit_string + '\n\n')

        # Write Phase definition to file
        element_names = ' '.join(trimmed_solution.element_names)
        element_names = element_names.replace('AR', 'Ar')
        species_names = trimmed_solution.species_names
        if species_qss:
            species_names = [spec for spec in species_names if spec not in species_qss]
            species_qss_names = build_species_string(species_qss, commented=True)
        species_names = build_species_string(species_names)

        if transport_model != 'Transport':
            transport_model_string = '     transport = \"$transport\", \n'
        else:
            transport_model_string = ''

        if kinetics:
            kinetics_string = '     kinetics = \"$kinetics\", \n'
        else:
            kinetics_string = ''

        if species_qss:
            phase_string = Template(
                'ideal_gas(name = \"$input_file_name_stripped\", \n' +
                '     elements = \"$elements\", \n' +
                '     species = """ $species""", \n' +
                '#    species_qss = """ $species_qss""", \n' +
                '     reactions = \"all\", \n' +
                # '     thermo = \"$thermo\", \n' + # not working for an obscure reason
                transport_model_string +
                kinetics_string +
                '     initial_state = state(temperature = $solution_temperature, '
                'pressure= $solution_pressure)   )       \n\n'
            )
        else:
            phase_string = Template(
                    'ideal_gas(name = \"$input_file_name_stripped\", \n' +
                    '     elements = \"$elements\", \n' +
                    '     species = """ $species""", \n' +
                    '     reactions = \"all\", \n' +
                    # '     thermo = \"$thermo\", \n' + # not working for an obscure reason
                    transport_model_string +
                    kinetics_string +
                    '     initial_state = state(temperature = $solution_temperature, '
                    'pressure= $solution_pressure)   )       \n\n'
            )

        if species_change_dict:
            # Name changing part
            species_names = species_names.split()
            for index, name in enumerate(species_names):
                if name in species_change_dict:
                    species_names[index] = species_change_dict[name]
            step = 10
            insert_index = 10
            while insert_index < index:
                species_names.insert(insert_index, '\n\t\t\t\t')
                insert_index += step + 1
            species_names = ' '.join(species_names)

        if species_qss:
            f.write(phase_string.substitute(
                elements=element_names,
                species=species_names,
                species_qss=species_qss_names,
                input_file_name_stripped=input_file_name_stripped,
                # thermo=thermo,
                transport=transport_model,
                kinetics=kinetics,
                solution_temperature=solution_temperature,
                solution_pressure=solution_pressure
            ))
        else:
            f.write(phase_string.substitute(
                    elements=element_names,
                    species=species_names,
                    input_file_name_stripped=input_file_name_stripped,
                    # thermo=thermo,
                    transport=transport_model,
                    kinetics=kinetics,
                    solution_temperature=solution_temperature,
                    solution_pressure=solution_pressure
            ))

        # Write species data to file
        section_break('Species data')
        for sp_index, name in enumerate(trimmed_solution.species_names):
            # joules/kelvin, boltzmann constant
            boltzmann = ct.boltzmann
            # 1 debye = d coulomb-meters
            debeye_conversion_constant = 3.33564e-30
            species = trimmed_solution.species(sp_index)
            name = str(trimmed_solution.species(sp_index).name)
            if name in species_change_dict:
                name = species_change_dict[name]
            nasa_coeffs = trimmed_solution.species(sp_index).thermo.coeffs
            replace_list_1 = {'{': '\"',
                              '}': '\"',
                              '\'': '',
                              ': ': ':',
                              '.0': "",
                              ',': '',
                              ' ': '  '}

            # build 7-coeff NASA polynomial array
            nasa_coeffs_1 = []
            for j, k in enumerate(nasa_coeffs):
                coeff = "{:.9e}".format(nasa_coeffs[j + 8])
                nasa_coeffs_1.append(coeff)
                if j == 6:
                    nasa_coeffs_1 = wrap_nasa(eliminate(str(nasa_coeffs_1),
                                                        {'\'': ""}))
                    break
            nasa_coeffs_2 = []
            for j, k in enumerate(nasa_coeffs):
                coeff = "{:.9e}".format(nasa_coeffs[j + 1])
                nasa_coeffs_2.append(coeff)
                if j == 6:
                    nasa_coeffs_2 = wrap_nasa(eliminate(
                        str(nasa_coeffs_2),
                        {'\'': ""}))
                    break
            # Species attributes from trimmed solution object
            composition = replace_multiple(
                str(species.composition),
                replace_list_1)

            nasa_range_1 = str([species.thermo.min_temp, nasa_coeffs[0]])
            nasa_range_2 = str([nasa_coeffs[0], species.thermo.max_temp])
            # check if species has defined transport data
            if bool(species.transport) is True:
                transport_geometry = species.transport.geometry
                diameter = str(species.transport.diameter * (10 ** 10))
                well_depth = str(species.transport.well_depth / boltzmann)
                polar = str(species.transport.polarizability * 10 ** 30)
                rot_relax = str(species.transport.rotational_relaxation)
                dipole = str(species.transport.dipole / debeye_conversion_constant)
                # create and fill string templates for each species
                if species.transport.dipole != 0:
                    species_string = Template(
                        'species(name = "$name",\n' +
                        '    atoms = $composition, \n' +
                        '    thermo = (\n' +
                        '       NASA(   $nasa_range_1, $nasa_coeffs_1  ),\n' +
                        '       NASA(   $nasa_range_2, $nasa_coeffs_2  )\n' +
                        '               ),\n'
                        '    transport = gas_transport(\n' +
                        '                   geom = \"$transport_geometry\",\n' +
                        '                   diam = $diameter, \n' +
                        '                   well_depth = $well_depth, \n' +
                        '                   polar = $polar, \n' +
                        '                   rot_relax = $rot_relax, \n' +
                        '                   dipole= $dipole) \n' +
                        '        )\n\n'
                    )
                    f.write(species_string.substitute(
                        name=name,
                        composition=composition,
                        nasa_range_1=nasa_range_1,
                        nasa_coeffs_1=nasa_coeffs_1,
                        nasa_range_2=nasa_range_2,
                        nasa_coeffs_2=nasa_coeffs_2,
                        transport_geometry=transport_geometry,
                        diameter=diameter,
                        well_depth=well_depth,
                        polar=polar,
                        rot_relax=rot_relax,
                        dipole=dipole
                    ))
                if species.transport.dipole == 0:
                    species_string = Template(
                        'species(name = "$name",\n'
                        '    atoms = $composition, \n'
                        '    thermo = (\n'
                        '       NASA(   $nasa_range_1, $nasa_coeffs_1  ),\n'
                        '       NASA(   $nasa_range_2, $nasa_coeffs_2  )\n'
                        '               ),\n'
                        '    transport = gas_transport(\n'
                        '                   geom = \"$transport_geometry\",\n'
                        '                   diam = $diameter, \n'
                        '                   well_depth = $well_depth, \n'
                        '                   polar = $polar, \n'
                        '                   rot_relax = $rot_relax) \n'
                        '        )\n\n'
                    )
                    f.write(species_string.substitute(
                        name=name,
                        composition=composition,
                        nasa_range_1=nasa_range_1,
                        nasa_coeffs_1=nasa_coeffs_1,
                        nasa_range_2=nasa_range_2,
                        nasa_coeffs_2=nasa_coeffs_2,
                        transport_geometry=transport_geometry,
                        diameter=diameter,
                        well_depth=well_depth,
                        polar=polar,
                        rot_relax=rot_relax,
                    ))
            if bool(species.transport) is False:
                species_string = Template(
                    'species(name = "$name",\n'
                    '    atoms = $composition, \n'
                    '    thermo = (\n'
                    '       NASA(   $nasa_range_1, $nasa_coeffs_1  ),\n'
                    '       NASA(   $nasa_range_2, $nasa_coeffs_2  )\n'
                    '               ),\n'
                    '        )\n\n'
                )
                f.write(species_string.substitute(
                    name=name,
                    composition=composition,
                    nasa_range_1=nasa_range_1,
                    nasa_coeffs_1=nasa_coeffs_1,
                    nasa_range_2=nasa_range_2,
                    nasa_coeffs_2=nasa_coeffs_2,
                ))

        # Write reactions to file

        section_break('Reaction Data')

        if species_qss and dummy:
            f.write("# Dummy reaction\n" + dummy + "\n\n")

        # write data for each reaction in the Solution Object
        for eq_index in range(len(trimmed_solution.reaction_equations())):
            equation_string = str(trimmed_solution.reaction_equation(eq_index))
            equation_object = trimmed_solution.reaction(eq_index)
            equation_type = type(equation_object).__name__
            m = str(eq_index + 1)
            if equation_type == 'ThreeBodyReaction':
                # trims efficiencies list
                efficiencies = equation_object.efficiencies
                trimmed_efficiencies = equation_object.efficiencies
                for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]
                arrhenius = build_arrhenius(equation_object, equation_type)
                replace_list_2 = {"{": "\"",
                                  "\'": "",
                                  ": ": ":",
                                  ",": " ",
                                  "}": "\""}

                efficiencies_string = replace_multiple(
                    str(trimmed_efficiencies),
                    replace_list_2)
                if species_qss:
                    template_string = '##  Reaction $m\n' \
                                      + '#three_body_reaction( \"$equation_string\",  $Arr'
                    if len(efficiencies) > 1:
                        template_string += ',\n#       efficiencies = $Efficiencies'
                    template_string += ') \n\n'

                    reaction_string = Template(template_string)
                else:
                    template_string = '#  Reaction $m\n' \
                                      + 'three_body_reaction( \"$equation_string\",  $Arr'
                    if len(efficiencies) > 1:
                        template_string += ',\n       efficiencies = $Efficiencies'
                    template_string += ') \n\n'

                    reaction_string = Template(template_string)
                if species_change_dict:
                    # Name changing part
                    efficiencies_string_2 = replace_multiple(efficiencies_string, {'"': ' " ', ':': ' : '})
                    efficiencies_string_temp = efficiencies_string_2.split()
                    for name in species_change_dict:
                        if name in efficiencies_string_temp:
                            for index, string in enumerate(efficiencies_string_temp):
                                if string == name:
                                    efficiencies_string_temp[index] = species_change_dict[name]
                    efficiencies_string = ' '.join(efficiencies_string_temp)

                    equation_string_temp = equation_string.split()
                    for name in species_change_dict:
                        if name in equation_string_temp:
                            for index, string in enumerate(equation_string_temp):
                                if string == name:
                                    equation_string_temp[index] = species_change_dict[name]
                    equation_string = ' '.join(equation_string_temp)

                f.write(reaction_string.substitute(
                    m=m,
                    equation_string=equation_string,
                    Arr=arrhenius,
                    Efficiencies=efficiencies_string
                ))
            if equation_type == 'ElementaryReaction':
                arrhenius = build_arrhenius(equation_object, equation_type)
                reaction_string = '#   Reaction $m\n' \
                                  'reaction( \"$equation_string\", $Arr'

                if equation_object.duplicate:
                    reaction_string += ',\n         options = \'duplicate\''
                    if equation_object.rate.pre_exponential_factor < 0:
                        reaction_string = reaction_string[:-1]
                        reaction_string += ', \'negative_A\']'

                if equation_object.orders:
                    order_string = ""
                    spec_keys = list(equation_object.orders.keys())
                    for spec in spec_keys:
                        order_string += f"{spec}:{equation_object.orders[spec]}"
                        if spec != spec_keys[-1]:
                            order_string += ' '

                    reaction_string += f',\n         order="{order_string}"'

                reaction_string += ')\n\n'

                if species_qss:
                    reaction_string = reaction_string.replace('#   Reaction', '##  Reaction')
                    reaction_string = reaction_string.replace('reaction', '#reaction')
                    reaction_string = reaction_string.replace('         options', '#        options')
                    reaction_string = reaction_string.replace('         order', '#        order')

                reaction_string = Template(reaction_string)

                if species_change_dict:
                    # Name changing part
                    equation_string_temp = equation_string.split()
                    for name in species_change_dict:
                        if name in equation_string_temp:
                            for index, string in enumerate(equation_string_temp):
                                if string == name:
                                    equation_string_temp[index] = species_change_dict[name]
                    equation_string = ' '.join(equation_string_temp)

                f.write(reaction_string.substitute(
                    m=m,
                    equation_string=equation_string,
                    Arr=arrhenius
                ))
            if equation_type == 'PlogReaction':
                # separate the different arrhenius
                rates = equation_object.rates
                coeff_sum = sum(equation_object.reactants.values())
                if species_qss:
                    reaction_string = '##  Reaction {0}\n#pdep_arrhenius( \"{1}\",\n'.format(m, equation_string)
                else:
                    reaction_string = '#  Reaction {0}\npdep_arrhenius( \"{1}\",\n'.format(m, equation_string)

                for index, values in enumerate(rates):
                    P = '{:.6E}'.format(rates[index][0])
                    if coeff_sum == 1:
                        A = str('{:.6E}'.format(rates[index][1].pre_exponential_factor))
                    if coeff_sum == 2:
                        A = str('{:.6E}'.format(rates[index][1].pre_exponential_factor * 10 ** 3))
                    if coeff_sum == 3:
                        A = str('{:.6E}'.format(rates[index][1].pre_exponential_factor * 10 ** 6))
                    b = '{:.6E}'.format(rates[index][1].temperature_exponent)
                    Ea = '{:.6E}'.format(rates[index][1].activation_energy / calories_constant)
                    end = ')\n\n' if values == rates[-1] else ',\n'
                    if species_qss:
                        arr_string = '#               [{0}, {1}, {2}, {3}]{4}'.format(P, A, b, Ea, end)
                    else:
                        arr_string = '               [{0}, {1}, {2}, {3}]{4}'.format(P, A, b, Ea, end)
                    reaction_string += arr_string
                if equation_object.duplicate:
                    reaction_string = reaction_string[:-3]
                    if species_qss:
                        reaction_string += ",\n#               options='duplicate')\n\n"
                    else:
                        reaction_string += ",\n               options='duplicate')\n\n"
                f.write(reaction_string)

            if equation_type == 'FalloffReaction':
                # trims efficiencies list
                efficiencies = equation_object.efficiencies
                trimmed_efficiencies = equation_object.efficiencies
                for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]

                kf = build_modified_arrhenius(equation_object, 'high')
                kf0 = build_modified_arrhenius(equation_object, 'low')
                replace_list_2 = {
                    "{": "\"",
                    "\'": "",
                    ": ": ":",
                    ",": " ",
                    "}": "\""
                }
                efficiencies_string = replace_multiple(
                    str(trimmed_efficiencies),
                    replace_list_2)
                if species_qss:
                    template_string = '##  Reaction $m\n' \
                                      + '#falloff_reaction( \"$equation_string\",\n' \
                                      + '#        kf = $kf,\n' \
                                      + '#        kf0   = $kf0'
                    if len(efficiencies) > 1:
                        template_string += ',\n#        efficiencies = $Efficiencies'

                else:
                    template_string = '#  Reaction $m\n' \
                                      + 'falloff_reaction( \"$equation_string\",\n' \
                                      + '        kf = $kf,\n' \
                                      + '        kf0   = $kf0'
                    if len(efficiencies) > 1:
                        template_string += ',\n        efficiencies = $Efficiencies'

                if equation_object.duplicate:
                    if species_qss:
                        template_string += ",\n#        options='duplicate'"
                    else:
                        template_string += ",\n        options='duplicate'"

                reaction_string = Template(template_string)

                if species_change_dict:
                    # Name changing part
                    efficiencies_string_2 = replace_multiple(efficiencies_string, {'"': ' " ', ':': ' : '})
                    efficiencies_string_temp = efficiencies_string_2.split()
                    for name in species_change_dict:
                        if name in efficiencies_string_temp:
                            for index, string in enumerate(efficiencies_string_temp):
                                if string == name:
                                    efficiencies_string_temp[index] = species_change_dict[name]
                    efficiencies_string = ' '.join(efficiencies_string_temp)

                    equation_string_temp = equation_string.split()
                    for name in species_change_dict:
                        if name in equation_string_temp:
                            for index, string in enumerate(equation_string_temp):
                                if string == name:
                                    equation_string_temp[index] = species_change_dict[name]
                    equation_string = ' '.join(equation_string_temp)

                f.write(reaction_string.substitute(
                    m=m,
                    equation_string=equation_string,
                    kf=kf,
                    kf0=kf0,
                    Efficiencies=efficiencies_string))

                # If optional Arrhenius data included:
                try:
                    falloff_str = build_falloff(equation_object)
                    f.write(falloff_str)
                except IndexError:
                    f.write(')\n\n')

    return output_file_name
