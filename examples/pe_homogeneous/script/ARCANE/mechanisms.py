"""Module mechanisms.py handles mechanism object, on which reduction operates"""

# Import statements
import os
import subprocess
import filecmp
import time

import cantera as ct
import ARCANE.ct2cti as ct2cti
import ARCANE.ct2ck as ct2ck
import ARCANE.custom.custom_kinetics as custom
import ARCANE.networks as networks
import ARCANE.kwdict as kwdict
import ARCANE.tools as tools
import ARCANE.display as display
import ARCANE.database as database

logger = display.Logger()
logger.set_log('logCompute')

kwdict = kwdict.Kwdict()


def init_mechanism_database(mechdir='mechs', overwrite=False):
    """Initialize mechanism folder
    For now, do nothing if it exists. Might include cleanup routines later on.

    Parameters
    ----------
    mechdir :
        name of folder that will contain all mechanisms used during reduction, default is 'mechs'
    overwrite :
        if True, folder is overwritten (Default value = False)


    Created: 17/11/13 [PP]
    Last modified: 18/10/02 [QC]

    """

    # Initialize case directory
    database.create_dir(mechdir, overwrite)

    # Set up mechdir for all instances of Mechanism
    Mechanism.mechdir = mechdir


def example_mechanism():
    """Generates an example Mechanism object based on GRI-Mech 3.0


    Returns
    -------
    example_mech : class :func:`~ARCANE.mechanisms.Mechanism` object

    """

    example_mech = Mechanism('gri30.cti', name='GRI-Mech3.0', kinetics=None)

    return example_mech


class Mechanism:

    # Generator
    def __init__(self, cti, name=None, parent=None, how=None, species_qss_names=[], skeletal=None,
                 kinetics=None, transport=None, f90=None, pdep=None, overwrite=False,
                 fit_order=0, start_temperature=300, end_temperature=2500, number_of_temperatures=6000,
                 fit_pressure=101325):
        """Creates and acts on mechanism objects used in reduction

        Parameters
        ----------
        cti :
            Cantera Solution object (can also be the cti name)
        name :
            user defined name of the mechanism
        parent :
            instance of mechanism that led to current instance
        how :
            how me was obtained from parent (rmS, rmR, QSS etc.)
        species_qss_names :
            list of QSS species names
        skeletal :
            skeletal mechanism for mechanisms with QSS
        kinetics :
            enables to specify 'custom' as kinetics
        transport :
            changes the transport model of the Solution object
        f90 :
            if specified, automatically compiles it
        pdep :
            if a pressure value is specified, converts the cti to one without pdep reactions
        overwrite :
            if True, the mechanism existing under this name will be replaced
        fit_order :
            order of the polynomial used for fitting the reverse arrhenius functions
        start_temperature :
            start temperature for the fit
        end_temperature :
            end temperature for the fit
        number_of_temperatures :
            number of temperatures for the fit
        fit_pressure :
            pressure used in the equilibrium constant computation


        Created: 17/11/13[PP]
        Last modified: 22/02/25 [QC]
        """

        if database.database_system != 'database':
            if not hasattr(self, 'mechdir'):
                init_mechanism_database()
                logger.debug("Initiating default mechanisms directory : mechs")

        # Cantera Solution object
        if type(cti) is str:
            if not cti.endswith('.cti') and not cti.endswith('.xml'):
                cti += '.cti'
            if not os.path.isfile(cti) and not cti == 'gri30.cti':
                if database.database_system != 'database':
                    path_to_cti = database.database_path + '/mechs/' + cti
                else:
                    path_to_cti = database.database_path + '/' + cti[:-4] + '/' + cti
                if os.path.isfile(path_to_cti):
                    cti = path_to_cti
                else:
                    logger.error('\nERROR! Specified cti file ' + cti + ' does not exist !')
                    quit()

            # Automatically gets the QSS species
            if not species_qss_names:
                try:
                    species_qss_names = tools.get_qss_species(cti)
                except Exception:
                    logger.debug('Specified mechanism present in a remote location')

            self.ctmech = ct.Solution(cti)
            self.reset()

            # Extracting the kinetics
            if cti != 'gri30.cti':
                kinetics = tools.get_kinetics(cti)

        else:
            self.ctmech = cti

        # # Resetting mechanism
        # self.reset()
        # Dynamic library status
        self.dynamic_lib_opened = False

        self.ctmech.TP = 300, ct.one_atm
        self.species_qss_names = species_qss_names

        # Checking if transport data is available for the species
        # If so, default transport will be set to 'Mix'
        transport_data_found = [True if spec.transport else False for spec in self.ctmech.species()]
        if all(transport_data_found):
            self.ctmech.transport_model = 'Mix'

        # Setting the transport model
        if transport:
            self.ctmech.transport_model = transport

        self.transport = self.ctmech.transport_model
        self.kinetics = kinetics
        if self.kinetics == 'custom' and not f90:
            f90 = 'auto'

        # Parameters for the backward coefficient fitting
        self.fit_order = fit_order
        self.start_temperature = start_temperature
        self.end_temperature = end_temperature
        self.number_of_temperatures = number_of_temperatures
        self.fit_pressure = fit_pressure

        # Key characteristics
        self.ns = self.ctmech.n_species
        self.species_names = self.ctmech.species_names
        self.species = self.ctmech.species()

        self.nr = self.ctmech.n_reactions
        self.reactions = self.ctmech.reactions()
        self.nqss = len(self.species_qss_names)

        self.elements_names = self.ctmech.element_names

        # Network data
        self.network = networks.Network(self)

        # Parent
        self.parent = parent

        # Skeletal
        self.skeletal = skeletal

        # Converting pdep reactions to elementary
        # Not functional
        if pdep:
            self.pdep_to_elementary(pressure=pdep)
            self.reactions = self.ctmech.reactions()

        # Determine if the input is the skeletal
        is_skeletal = True
        for spec in self.species_qss_names:
            if not spec in self.species_names:
                is_skeletal = False
                break

        # Skeletal mechanism reconstruction
        if not self.skeletal:
            if self.nqss > 0:
                if is_skeletal:
                    self.skeletal = Mechanism(cti, parent=parent, how=how, skeletal=skeletal, transport=transport)
                    self.species_names = [spec for spec in self.species_names if spec not in self.species_qss_names]
                    self.ns = self.ns - self.nqss
                elif type(cti) == str:
                    tools.skeletal_from_custom_cti(cti, 'retrieved_skeletal.cti')

                    # Avoiding bug looking for a mixture_database.dat file for skeletal
                    if transport == 'AVBP':
                        transport_skeletal = None
                    else:
                        transport_skeletal = transport

                    self.skeletal = Mechanism('retrieved_skeletal.cti', transport=transport_skeletal, how='retrieved')
                    subprocess.call('rm -f retrieved_skeletal.cti', shell=True)
            else:
                self.skeletal = self

        # Lists of Species and Reactions objects as attributes
        self.species = [spec for spec in self.ctmech.species() if spec.name in self.species_names]

        # Set child of parent as self
        # parent.child = self

        # Child
        # TODO - Handle Mechanism instance having several children
        # self.child = None
        # Initialize array to keep track of multiple children.
        # self.childrenlist = ()

        # How it was obtained
        self.how = how

        # Timescales data
        self.timescales = {}

        # Construct name of mechanism.
        index = 0
        if self.nqss > 0:

            if self.skeletal:
                self.nr = self.skeletal.nr
                self.reactions = self.skeletal.ctmech.reactions()
                self.species_qss = [spec for spec in self.skeletal.species if spec.name in self.species_qss_names]

            elif type(cti) == str:
                self.nr = tools.get_number_of_reactions(cti)
                self.reactions = []

            else:
                self.nr = 1
                self.reactions = []

            if not parent:
                self.parent = self.skeletal

            prename = "S" + str(self.ns) + "R" + str(self.nr) \
                      + "QSS" + str(self.nqss)
            self.species_qss_names = species_qss_names
        else:
            prename = "S" + str(self.ns) + "R" + str(self.nr)
            self.reactions = self.ctmech.reactions()

        # Number of reversible reactions
        if self.species_qss_names and self.skeletal:
            self.nr_reverse = sum([1 for i in range(self.nr) if self.skeletal.ctmech.is_reversible(i)])
        else:
            self.nr_reverse = sum([1 for i in range(self.nr) if self.ctmech.is_reversible(i)])

        if not name:
            input_name = False
            # Iteratively check that it does not exists already
            name = prename + "_" + str(index)
            if database.database_system != 'database':
                while os.path.isfile(self.mechdir + "/" + name + '.cti'):
                    index += 1
                    name = prename + "_" + str(index)
            else:
                while os.path.isdir(database.database_path + "/" + name):
                    index += 1
                    name = prename + "_" + str(index)
        else:
            input_name = True
            if name == 'auto':
                name = cti.replace('.cti', '')

        self.name = name

        if database.database_system == 'database':
            # Creates mechanism directory if no already there (do not overwrite)
            self.mechdir = database.database_path + "/" + name
            database.create_dir(self.mechdir, overwrite=overwrite)

        # Print it to file
        self.path = self.mechdir + "/" + str(self.name) + ".cti"
        if self.nqss > 0:
            if self.skeletal:
                custom.cantera_qss_solution(self.skeletal.ctmech, species_qss_names,
                                            qss_cti_file=self.path, transport=transport,
                                            output_solution=False)
            else:
                custom.cantera_qss_solution(self.ctmech, species_qss_names,
                                            qss_cti_file=self.path, transport=transport,
                                            output_solution=False)
        else:
            ct2cti.write(self.ctmech, self.path, kinetics=self.kinetics, transport=transport)

        if not input_name:
            # If the mechanism are identical, do not copy
            new_file = self.mechdir + "/" + str(self.name) + ".cti"
            for step in range(index):
                name = prename + "_" + str(step)
                if database.database_system == 'database':
                    original_file = database.database_path + "/" + name + "/" + str(name) + ".cti"
                else:
                    original_file = self.mechdir + "/" + str(name) + ".cti"
                if filecmp.cmp(original_file, new_file):
                    subprocess.call('rm -f ' + new_file, shell=True)
                    if database.database_system == 'database':
                        subprocess.call('rm -rf ' + self.mechdir, shell=True)
                        self.mechdir = database.database_path + "/" + name
                    self.name = name
                    self.path = self.mechdir + "/" + str(self.name) + ".cti"
                    break

        # Compiles the f90 file
        if f90:
            self.kinetics = 'custom'
            if f90 != 'force_write':
                if not os.path.isfile(f90) or f90 == 'auto':
                    if os.path.isfile(self.mechdir + "/" + str(self.name) + ".f90"):
                        f90 = self.mechdir + "/" + str(self.name) + ".f90"
                    else:
                        if f90 == 'write' or f90 == 'auto':
                            f90 = self.mechdir + "/" + str(self.name) + ".f90"
                            custom.print_fortran(self, f90, self.species_qss_names,
                                                 fit_order=self.fit_order,
                                                 start_temperature=self.start_temperature,
                                                 end_temperature=self.end_temperature,
                                                 number_of_temperatures=self.number_of_temperatures,
                                                 pressure=self.fit_pressure)

                            logger.info('Fortran file automatically generated: ' + f90)
                        else:
                            logger.error('\nERROR! Specified f90 file ' + f90 + ' does not exists !')
                            quit()
                elif f90 != self.mechdir + "/" + str(self.name) + ".f90" and input_name:
                    subprocess.call('cp ' + f90 + ' ' + self.mechdir + "/" + str(self.name) + ".f90", shell=True)
                    f90 = self.mechdir + "/" + str(self.name) + ".f90"
            else:
                f90 = self.mechdir + "/" + str(self.name) + ".f90"

                custom.print_fortran(self, f90, self.species_qss_names,
                                     fit_order=self.fit_order,
                                     start_temperature=self.start_temperature,
                                     end_temperature=self.end_temperature,
                                     number_of_temperatures=self.number_of_temperatures,
                                     pressure=self.fit_pressure)

                logger.info('Fortran file automatically generated: ' + f90)

            self.f90 = f90

            filesize = os.path.getsize(self.f90)
            if filesize == 0:
                logger.warning(self.f90 + ' file was empty and will be rewritten')
                self.overwrite_fortran()

            elif not self.check_f90():
                logger.warning(self.f90 + ' file was corrupted and will be rewritten')
                self.overwrite_fortran()

        else:
            self.f90 = None

        # Nickname for legend in plotting
        if type(cti) is str and not input_name:
            self.nickname = cti.split('/')[-1]
            self.nickname = self.nickname.split('.')[0]
        else:
            self.nickname = self.name

        self.colour = None
        self.marker = None
        self.line_style = None

        # File name
        self.path = os.path.abspath(self.mechdir + "/" + self.name + '.cti')

        # Inert species
        self.inert = []

        # Species keyword dictionaries for post-processing
        kwdict_names_extension = {}
        kwdict_units_extension = {}
        for spec in self.species_names:
            # Mass fraction
            kwdict_names_extension[spec] = []
            for suffix in kwdict.names['Y']:
                kwdict_names_extension[spec].append(spec + ' ' + suffix)
            kwdict_names_extension[spec].append(spec)
            kwdict_names_extension['Y_' + spec] = kwdict_names_extension[spec]

            # Mole fraction
            kwdict_names_extension['X_' + spec] = []
            for suffix in kwdict.names['X']:
                kwdict_names_extension['X_' + spec].append(spec + ' ' + suffix)
            kwdict_names_extension['X_' + spec].append('X_' + spec)

            # Concentration
            kwdict_names_extension['c_' + spec] = []
            for suffix in kwdict.names['c']:
                kwdict_names_extension['c_' + spec].append(spec + ' ' + suffix)
            kwdict_names_extension['c_' + spec].append('c_' + spec)

            # Production rates
            kwdict_names_extension['cdot_' + spec] = []
            for suffix in kwdict.names['cdot']:
                kwdict_names_extension['cdot_' + spec].append(spec + ' ' + suffix)
            kwdict_names_extension['cdot_' + spec].append('cdot_' + spec)

            kwdict_units_extension[kwdict_names_extension[spec][0]] = kwdict.units['Mass fraction']
            kwdict_units_extension[kwdict_names_extension['Y_' + spec][0]] = kwdict.units['Mass fraction']
            kwdict_units_extension[kwdict_names_extension['X_' + spec][0]] = kwdict.units['Mole fraction']
            kwdict_units_extension[kwdict_names_extension['c_' + spec][0]] = kwdict.units['Concentration']
            kwdict_units_extension[kwdict_names_extension['cdot_' + spec][0]] = kwdict.units['Production Rates']

        if self.f90 and self.skeletal:
            reaction_to_iterate = self.skeletal.ctmech.reaction_equations()
        else:
            reaction_to_iterate = self.ctmech.reaction_equations()

        for id_reac, reac in enumerate(reaction_to_iterate):
            # Reactions rates
            kwdict_names_extension['w_' + str(id_reac)] = []
            for suffix in kwdict.names['w']:
                kwdict_names_extension['w_' + str(id_reac)].append(str(id_reac) + ' ' + suffix)
            kwdict_names_extension['w_' + str(id_reac)].append('w_' + str(id_reac))
            kwdict_units_extension[kwdict_names_extension['w_' + str(id_reac)][0]] = kwdict.units['Reaction Rates']

        self.kwdict_names_extension = kwdict_names_extension
        self.kwdict_units_extension = kwdict_units_extension

        # Updating kwdict with the species_names
        kwdict.update('names', self.kwdict_names_extension)
        kwdict.update('units', self.kwdict_units_extension)

        # Parameters of the merging
        self.base_mechanism = self
        self.additional_mechanisms = []
        self.n_add_mechs = 0
        self.merged_mechanism = None

    def has_qss(self):
        """Finds if the mechanism has QSS species

        Returns
        -------
        has_qss : bool
            true if the mechanism has qss
        skeletal_mechanism : class `ARCANE.Mechanism` object
            skeletal mecanism, self if has not qss
        species_qss_names : list
            qss species in a list

        """
        if self.species_qss_names:
            has_qss = True
            skeletal_mechanism = self.skeletal

            # QSS instantiation
            species_qss_names = self.species_qss_names
        else:
            has_qss = False
            skeletal_mechanism = self

            # QSS instantiation
            species_qss_names = []

        return has_qss, skeletal_mechanism, species_qss_names

    def rename(self, new_name):
        """Changes the name of the mechanism

        Parameters
        ----------
        new_name :
            new name for the mechanism

        """
        if os.path.isdir(database.database_path + "/" + new_name):
            logger.error(f'ERROR! New name {new_name} already exists, you should choose another one')
            return

        if database.database_system != 'database':
            logger.warning('WARNING! Renaming for classical storage system is not coded yet')
            return

        old_name = self.name.replace('.cti', '')

        # New mechanism directory name
        new_mech_dir = database.database_path + "/" + new_name

        # Renaming the old directory
        subprocess.call('mv ' + self.mechdir + ' ' + new_mech_dir, shell=True)

        # New .cti name
        old_cti_name = old_name + '.cti'
        new_cti_name = new_name + '.cti'

        # Renaming the .cti
        subprocess.call('mv ' + new_mech_dir + '/' + old_cti_name + ' '
                        + new_mech_dir + '/' + new_cti_name, shell=True)

        # Changing the attributes values
        self.mechdir = new_mech_dir
        self.name = new_name
        self.nickname = new_name
        self.path = self.mechdir + '/' + new_cti_name

        if self.f90:
            # New f90 name
            old_f90_name = old_name + '.f90'
            new_f90_name = new_name + '.f90'

            # Renaming the f90
            subprocess.call('mv ' + new_mech_dir + '/' + old_f90_name + ' '
                            + new_mech_dir + '/' + new_f90_name, shell=True)
            self.f90 = self.mechdir + '/' + new_f90_name

    def overwrite_fortran(self):
        """Rewrites the f90 file with the current parameters

        """

        custom.print_fortran(self, self.f90, self.species_qss_names,
                             fit_order=self.fit_order,
                             start_temperature=self.start_temperature,
                             end_temperature=self.end_temperature,
                             number_of_temperatures=self.number_of_temperatures,
                             pressure=self.fit_pressure)

    def compile(self):
        """Compiles the mechanism if this is an ARC

        """

        if self.f90:
            logger.debug('Compiling ' + self.f90)
            custom.compile_fortran(self.f90)
            self.ctmech = ct.Solution(self.path)
            self.dynamic_lib_opened = True
        else:
            logger.debug('Mechanism object has not f90 attribute, nothing will be done')

    def reset(self):
        """Resets the dynamic library

        """

        if hasattr(self, 'name'):
            logger.debug(f'Resetting dynamic library of {self.name}')
        if hasattr(self.ctmech, 'reset_custom'):
            self.ctmech.reset_custom()
            self.dynamic_lib_opened = False
            time.sleep(1)

    def check_f90(self, debug=False):
        """Checks the f90 file and throws an error if it is not valid


        Returns
        -------
        exit_status : bool
            status, True if file is OK, False else

        """

        #  Test to know if subroutine has the right name to do the compiling
        line_custom = None
        exit_status = True

        f = open(self.f90, 'r')
        for lines in f.readlines():
            if 'subroutine' in lines.lower() and '(p,t,y,wdot)' in lines.replace(" ", "").lower():
                line_custom = lines
                break
        f.close()

        if not line_custom:
            logger.debug("The main subroutine of the f90 file should be written : "
                         "'subroutine customkinetics (P, T, y, wdot)' to work")
            exit_status = False
        else:
            if 'customkinetics' not in line_custom:
                logger.debug("The main subroutine of the f90 file should be written : "
                             "'subroutine customkinetics (P, T, y, wdot)' to work")
                exit_status = False

        return exit_status

    def copy(self, new_name=None, directory=None):
        """Copy the Mechanism object

        Parameters
        ----------
        new_name :
            name of the new mechanism (Default value = None)
        directory :
            directory in which the file will be stored (Default value = None)

        Returns
        -------
        duplicate_mechanism : class `ARCANE.mechanisms.Mechanism` object
            Copy of a class `ARCANE.mechanisms.Mechanism` object

        """
        if self.transport == 'Transport':
            transport_to_duplicate = None
        else:
            transport_to_duplicate = self.transport

        mech_dir_save = self.mechdir
        if directory:
            self.mechdir = directory

        duplicate_mechanism = Mechanism(self.path, name=new_name, parent=self.parent,
                                        how='copy', species_qss_names=self.species_qss_names,
                                        skeletal=self.skeletal, kinetics=self.kinetics,
                                        transport=transport_to_duplicate, f90=self.f90)

        self.mechdir = mech_dir_save

        if directory:
            subprocess.call('cp ' + duplicate_mechanism.path + ' ' + directory + '/', shell=True)
            if self.f90:
                subprocess.call('cp ' + duplicate_mechanism.f90 + ' ' + directory + '/', shell=True)

        return duplicate_mechanism

    def change_species_name(self, species_to_change=None, new_name=None, characters_change_dict=None,
                            file_name='modified_cti.cti'):
        """
        Changes the name of the species in the mechanism

        Parameters
        ----------
        species_to_change : list
            list of the species to be changed, if nothing is specified, the change will be made on all names
        new_name : list
            list of the new species names, should match the length of the previous argument
        characters_change_dict : dict
            dictionary providing the correspondence of characters to be modified in the list species_to_change
        file_name : str
            name of the new cti file

        Returns
        -------

        """
        # Checking argument type
        if species_to_change:
            if type(species_to_change) != list:
                species_to_change = list(species_to_change)
        else:
            species_to_change = self.species_names

        all_species_objects = self.species

        loaded_cti = open(self.path)
        string_cti = loaded_cti.read()

        if new_name:
            if type(new_name) != list:
                new_name = list(new_name)
            if len(new_name) != len(species_to_change):
                logger.error('The number of new names must correspond to the number of species to change. '
                             'Nothing will be done.')
                return

            for index, spec in enumerate(species_to_change):
                string_cti = string_cti.replace(spec, new_name(index))

        if characters_change_dict:
            modified_species = []
            for spec in species_to_change:
                new_spec = spec
                for key in characters_change_dict:
                    if key in spec:
                        new_spec = new_spec.replace(key, characters_change_dict[key])

                # Change in cti
                string_cti = string_cti.replace(spec, new_spec)

        test_cti = open(file_name, 'w+')
        test_cti.write(string_cti)

        loaded_cti.close()
        test_cti.close()

    def remove_species(self, species_to_discard_list, name=None, auto_remove=False):
        """Removing species from the mechanism

        Parameters
        ----------
        species_to_discard_list :
            list of species to discard
        name :
            name of the new Mechanism object (Default value = None)
        auto_remove :
            if True, removes species that are no longer used in a reaction (Default value = False)

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object with discarded species

        """

        # Converting potential string input to list
        if not isinstance(species_to_discard_list, list):
            species_to_discard_list = list(species_to_discard_list)

        if self.f90:

            phase = cantera_reduced_solution(self.skeletal.ctmech,
                                             species_to_discard=species_to_discard_list,
                                             diluent=self.inert,
                                             auto_remove=auto_remove,
                                             transport=self.transport)

            new_species_qss_names = [spec for spec in self.species_qss_names if spec not in species_to_discard_list]

            new_mechanism = Mechanism(phase, species_qss_names=new_species_qss_names,
                                      kinetics='custom', f90='write', name=name, parent=self)

        else:

            phase = cantera_reduced_solution(self.ctmech,
                                             species_to_discard=species_to_discard_list,
                                             diluent=self.inert,
                                             auto_remove=auto_remove,
                                             transport=self.transport)

            new_mechanism = Mechanism(phase, how='species removal', name=name, parent=self)
            new_mechanism.inert = self.inert

        return new_mechanism

    def add_species(self, species_to_add_list, name=None, auto_remove=False):
        """Adds species to the mechanism

        Parameters
        ----------
        species_to_add_list :
            list of species objects to discard
        name :
            name of the new Mechanism object (Default value = None)
        auto_remove :
            if True, removes species that are no longer used in a reaction (Default value = False)

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object with added species

        """

        # Converting potential string input to list
        if not isinstance(species_to_add_list, list):
            species_to_add_list = [species_to_add_list]

        if self.f90:

            phase = cantera_reduced_solution(self.skeletal.ctmech,
                                             species_to_add=species_to_add_list,
                                             diluent=self.inert,
                                             auto_remove=auto_remove,
                                             transport=self.transport)

            new_species_qss_names = [spec for spec in self.species_qss_names if spec not in species_to_add_list]

            new_mechanism = Mechanism(phase, species_qss_names=new_species_qss_names,
                                      kinetics='custom', f90='write', name=name, parent=self)

        else:

            phase = cantera_reduced_solution(self.ctmech,
                                             species_to_add=species_to_add_list,
                                             diluent=self.inert,
                                             auto_remove=auto_remove,
                                             transport=self.transport)

            new_mechanism = Mechanism(phase, how='species addition', name=name, parent=self)
            new_mechanism.inert = self.inert

        return new_mechanism

    def remove_reactions(self, reactions_list, name=None, auto_remove=False):
        """Removing reactions from the mechanism

        Parameters
        ----------
        reactions_list :
            list of reactions to discard (can be indices)
        name :
            name of the new Mechanism object (Default value = None)
        auto_remove :
            if True, removes species that are no longer used in a reaction (Default value = False)

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object with discarded reactions

        """

        # Converting potential string input to list
        if not isinstance(reactions_list, list):
            reactions_list = list(reactions_list)

        # Converting reactions indices to reaction objects
        for index, reac in enumerate(reactions_list):
            if isinstance(reac, int):
                reactions_list[index] = self.reactions[index]

        if self.f90:
            phase = cantera_reduced_solution(self.skeletal.ctmech,
                                             reactions_to_discard=reactions_list,
                                             diluent=self.inert,
                                             auto_remove=auto_remove,
                                             transport=self.transport)

            new_mechanism = Mechanism(phase, species_qss_names=self.species_qss_names,
                                      kinetics='custom', f90='write', name=name, parent=self)

        else:

            phase = cantera_reduced_solution(self.ctmech,
                                             reactions_to_discard=reactions_list,
                                             diluent=self.inert,
                                             auto_remove=auto_remove,
                                             transport=self.transport)

            new_mechanism = Mechanism(phase, how='reactions removal', name=name, parent=self)
            new_mechanism.inert = self.inert

        return new_mechanism

    def add_reaction(self, reactions_list, name=None):
        """Adding reactions to the mechanism

        Parameters
        ----------
        reactions_list :
            list of reactions to discard (can be indices)
        name :
            name of the new Mechanism object (Default value = None)

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object with added reactions

        """

        # Converting potential string input to list
        if not isinstance(reactions_list, list):
            reactions_list = list(reactions_list)

        if self.f90:
            phase = cantera_reduced_solution(self.skeletal.ctmech,
                                             reactions_to_add=reactions_list,
                                             diluent=self.inert,
                                             transport=self.transport)

            new_mechanism = Mechanism(phase, species_qss_names=self.species_qss_names,
                                      kinetics='custom', f90='write', name=name, parent=self)

        else:

            phase = cantera_reduced_solution(self.ctmech,
                                             reactions_to_add=reactions_list,
                                             diluent=self.inert,
                                             transport=self.transport)

            new_mechanism = Mechanism(phase, how='reactions addition', name=name, parent=self)
            new_mechanism.inert = self.inert

        return new_mechanism

    def add_qss(self, species_list, name=None):
        """Adds a list of QSS species to the mechanism

        Parameters
        ----------
        species_list :
            list of species to be added to the QSS species
        name :
            name of the new Mechanism object (Default value = None)

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object with new QSS list

        """

        species_qss_names = []
        if not isinstance(species_list, list):
            species_list = [species_list]

        if self.nqss == 0:
            species_qss_names = species_list
            ctmech = self.ctmech
        else:
            species_qss_names += species_list
            # Avoiding doubles
            species_qss_names = list(set(species_qss_names))
            ctmech = self.skeletal.ctmech

        new_mechanism = Mechanism(ctmech, species_qss_names=species_qss_names,
                                  kinetics='custom', f90='write', name=name, parent=self)

        return new_mechanism

    def add_mechanism(self, additional_mechanisms):
        """Adds a mechanism to the existing one

        Parameters
        ----------
        additional_mechanisms :
            list of additional Mechanism objects


        """

        if self.nqss > 0:
            logger.warning('\n!!! WARNING !!! You are trying to merge a mechanism into an ARC.')
            logger.warning('Only transported species will be merged and the resulting cti will be useless')
            logger.warning('Retrieve the skeletal mechanism (self.skeletal) if you want a suitable mechanism')

        if not isinstance(additional_mechanisms, list):
            additional_mechanisms = [additional_mechanisms]

        self.additional_mechanisms += additional_mechanisms
        self.n_add_mechs += 1

    def cut_mechanism(self, criteria={}, weights=None, force_keep=[], force_remove=[],
                      auto_remove=True, name=None, cti_name=None):
        """Remove each species that has more atoms than criteria[atom]

        Parameters
        ----------
        criteria :
            dictionary of minimum number of atoms (Default value = {})
        weights :
            list of molecular weights to target (Default value = None)
        force_keep :
            species that will be kept even if meeting the criteria (Default value = [])
        force_remove :
            species that will be removed even if not meeting the criteria (Default value = [])
        auto_remove :
            if True, removes species that are no longer used in a reaction (Default value = True)
        name :
            name of the new Mechanism object (Default value = None)
        cti_file :
            name of the ouput cti
        cti_name :
             (Default value = None)

        Returns
        -------
        new_mech : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object following the given criteria

        """
        ctmech = self.ctmech

        species = ctmech.species_names
        species_to_discard = force_remove

        for spec in species:
            for atom in criteria:
                if criteria[atom] >= 0:
                    if ctmech.n_atoms(spec, atom) > criteria[atom]:
                        species_to_discard.append(spec)
                else:
                    if ctmech.n_atoms(spec, atom) < - criteria[atom]:
                        species_to_discard.append(spec)

            if weights:
                for weight in weights:
                    if weight >= 0:
                        if ctmech.molecular_weights[species.index(spec)] > weight:
                            species_to_discard.append(spec)
                    else:
                        if ctmech.molecular_weights[species.index(spec)] < - weight:
                            species_to_discard.append(spec)

        species_to_discard = [spec for spec in species_to_discard if spec not in force_keep]
        species_to_discard = list(dict.fromkeys(species_to_discard))
        new_ctmech = cantera_reduced_solution(ctmech=ctmech,
                                              species_to_discard=species_to_discard,
                                              diluent=force_keep,
                                              reduced_cti_file=cti_name,
                                              auto_remove=auto_remove)

        new_mech = Mechanism(new_ctmech, how='cut', name=name)

        head, tail = display.head_tail_charac(text='Removed species')
        logger.info(head)
        logger.info('\n'.join(species_to_discard))
        logger.info('\nPrevious mechanism : ns = ' + str(self.ns) + ', nr = ' + str(self.nr))
        logger.info('New mechanism : ns = ' + str(new_mech.ns) + ', nr = ' + str(new_mech.nr))
        if cti_name:
            logger.info('New mechanism saved as : ' + cti_name)

        return new_mech

    def extract_mechanism(self, target_species_list, name=None):
        """Extracts the sub-mechanism involving the targeted species

        Parameters
        ----------
        target_species_list :
            list of species to be retrieved by the sub-mechanism
        name :
            name of the new Mechanism object (Default value = None)

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object with the targeted species

        """

        if not isinstance(target_species_list, list):
            target_species_list = [target_species_list]

        head, tail = display.head_tail_charac('Extracting mechanism')
        logger.info(head)

        species_involved = set()
        reactions_to_keep = []

        species = self.ctmech.species()
        species_names = self.species_names
        reactions = self.ctmech.reactions()

        for spec in target_species_list:
            if spec not in species_names:
                logger.error('\nSpecies ' + spec + ' is not in the mechanism, it will be ignored')
                target_species_list.remove(spec)

        for reac in reactions:
            for spec in target_species_list:
                if spec in reac.reactants or spec in reac.products:
                    species_involved.update(reac.reactants)
                    species_involved.update(reac.products)

                    reactions_to_keep.append(reac)
                    break

        species_involved = list(species_involved)
        species_to_keep = [spec for spec in species if spec.name in species_involved]

        logger.info('\nThe following species have been extracted: ' + ', '.join(species_involved))
        logger.info('\nThe following reactions have been extracted:')
        for item in reactions_to_keep:
            logger.info(item)

        phase = ct.Solution(thermo='IdealGas',
                            transport=self.transport,
                            kinetics='GasKinetics',
                            species=species_to_keep,
                            reactions=reactions_to_keep)

        new_mechanism = Mechanism(phase, how='extraction', name=name)

        logger.info(tail)

        return new_mechanism

    def merge_species(self, name_flags=[], keep_isomers=False):
        """Merges the species objects of two mechanisms

        Parameters
        ----------
        name_flags :
            list of characters chains flagging isomer species to keep ('*' for example) (Default value = [])
        keep_isomers :
            if True, identical species will be kept as is (Default value = False)

        Returns
        -------
        merged_species : list
            list of merged species
        changed_names_dict : dict
            dictionary of names correspondence

        """

        if self.additional_mechanisms == 0:
            logger.error('There is no supplementary mechanism to merge into the base one')
            quit()

        base_ctmech = self.base_mechanism.ctmech

        base_species = base_ctmech.species()
        merged_species = base_ctmech.species()
        changed_names_dict = {}

        for index, addi_mech in enumerate(self.additional_mechanisms):

            logger.info('\n==> Merging species from ' + addi_mech.name)

            new_mechanism = addi_mech
            new_ctmech = new_mechanism.ctmech
            new_species = new_ctmech.species()

            added_species = []
            changed_names_dict_local = {}

            for species in new_species:
                add_species = True
                species_to_add = species
                for b_species in base_species:
                    if tools.same_species(species, b_species):
                        if species.name == b_species.name:
                            add_species = False
                            break
                        else:
                            if not keep_isomers:
                                keep_isomer = False
                                if name_flags:
                                    for flag in name_flags:
                                        if flag in species.name:
                                            keep_isomer = True
                                            continue

                                if not keep_isomer:
                                    add_species = False
                                    changed_names_dict_local[species.name] = b_species.name
                    else:
                        if species.name == b_species.name:
                            if keep_isomers:
                                new_name = species.name + 'ISOMER'
                                species_to_add = tools.create_species(new_name, species)
                                changed_names_dict_local[species.name] = new_name
                            else:
                                add_species = False

                if add_species:
                    added_species.append(species_to_add)

            merged_species += added_species
            base_species = merged_species

            if changed_names_dict_local:
                logger.info('\nThe following species will be renamed to match the base mechanism:')
                for spec in changed_names_dict_local:
                    logger.info(spec + ' --> ' + changed_names_dict_local[spec])

            list_of_names = [spec.name for spec in added_species]
            if list_of_names:
                logger.info('\nThe following species will be added: ' + ', '.join(list_of_names))
            else:
                logger.info('\nNothing to be added from this mechanism')

            changed_names_dict.update(changed_names_dict_local)

        return merged_species, changed_names_dict

    def merge_reactions(self, changed_names_dict=None):
        """Merges the reactions objects of two mechanisms

        Parameters
        ----------
        changed_names_dict :
            dictionary of species names that need to be changed (Default value = None)

        Returns
        -------
        merged_reactions : list
            list of merged reactions

        """

        if self.additional_mechanisms == 0:
            logger.error('There is no supplementary mechanism to merge into the base one')
            quit()

        base_ctmech = self.base_mechanism.ctmech

        base_reactions = base_ctmech.reactions()
        merged_reactions = base_ctmech.reactions()

        for index, addi_mech in enumerate(self.additional_mechanisms):

            added_reactions = []

            logger.info('\n==> Merging reactions from ' + addi_mech.name)

            new_mechanism = addi_mech
            new_ctmech = new_mechanism.ctmech
            new_reactions = new_ctmech.reactions()

            if changed_names_dict:
                new_reactions = tools.change_reaction_by_species(new_ctmech, list(changed_names_dict.keys()),
                                                                 new_names_list=list(changed_names_dict.values()))

            for new_reac in new_reactions:
                add_reac = True

                for base_reac in base_reactions:
                    if tools.same_reaction(new_reac, base_reac):
                        add_reac = False
                        continue

                if add_reac:
                    added_reactions.append(new_reac)

            merged_reactions += added_reactions
            base_reactions = merged_reactions

            if added_reactions:
                logger.info('\nThe following reactions will be added to the base mechanism:')
                for item in added_reactions:
                    logger.info(item.equation)
            else:
                logger.info('\nNothing to be added from this mechanism')

        return merged_reactions

    def merge_mechanisms(self, name=None, keep_isomers=False, name_flags=[]):
        """Merges 2 mechanisms

        Parameters
        ----------
        name :
            name of the new Mechanism object (Default value = None)
        keep_isomers :
            if True, identical species will be kept as is (Default value = False)

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object merging two mechanisms

        """

        head, tail = display.head_tail_charac('Merging mechanisms')
        logger.info(head)

        merged_species, changed_names_dict = self.merge_species(keep_isomers=keep_isomers, name_flags=name_flags)
        merged_reactions = self.merge_reactions(changed_names_dict=changed_names_dict)

        logger.info(tail)

        new_ctmech = ct.Solution(thermo='IdealGas',
                                 transport=self.transport,
                                 kinetics='GasKinetics',
                                 species=merged_species,
                                 reactions=merged_reactions)

        new_mechanism = Mechanism(new_ctmech, how='merging', name=name)

        logger.info('Merged mechanism : ns = ' + str(new_mechanism.ns) + ', nr = ' + str(new_mechanism.nr))

        return new_mechanism

    def pdep_to_elementary(self, name=None, pressure=1e5):
        """Converting a Plog reaction to an elementary reaction

        Parameters
        ----------
        name :
            name of the new Mechanism object (Default value = None)
        pressure :
            pressure for which the closest Arrhenius will be taken

        Returns
        -------
        new_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            new class :func:`~ARCANE.mechanisms.Mechanism` object without Plog reactions (Default value = 1e5)

        """

        ctmech = self.ctmech

        reactions = ctmech.reactions()

        new_reactions = []

        for reaction in reactions:

            # Extracting data from the old reaction
            reactants = reaction.reactants
            products = reaction.products

            reaction_type = reaction.reaction_type

            if reaction_type == 5:

                # New pressure dependent reaction
                new_reaction = ct.ElementaryReaction(reactants, products)

                pressure_list = []

                for rate in reaction.rates:
                    pressure_list.append(rate[0])

                index_bar = pressure_list.index(min(pressure_list, key=lambda x: abs(x - pressure)))

                new_reaction.rate = reaction.rates[index_bar][1]

                new_reaction.ID = reaction.ID
                new_reaction.duplicate = reaction.duplicate
                new_reaction.orders = reaction.orders
                new_reaction.reversible = reaction.reversible

            else:
                new_reaction = reaction

            new_reactions.append(new_reaction)

        # Creates new Solution object
        new_ctmech = ct.Solution(thermo='IdealGas',
                                 transport=self.transport,
                                 kinetics='GasKinetics',
                                 species=self.species,
                                 reactions=new_reactions)

        new_mechanism = Mechanism(new_ctmech, name=name, parent=self.parent,
                                  how='pdep suppression')

        return new_mechanism

    def write(self, file_name=None):
        """Writes the mechanism in cti format

        Parameters
        ----------
        file_name :
            name of the output mechanism (Default value = None)

        """

        if not file_name:
            file_name = self.name + '.cti'

        if self.nqss and self.skeletal:
            custom.cantera_qss_solution(self.skeletal.ctmech, self.species_qss_names,
                                        qss_cti_file=file_name, transport=self.transport,
                                        output_solution=False)
            logger.info(f".cti file written as {file_name}")
        else:
            ct2cti.write(self.ctmech, file_name, species_qss=self.species_qss_names)
            logger.info(f".cti file written as {file_name}")

        if self.f90:
            if file_name == self.name + '.cti':
                file_name = file_name.replace('.cti', '.f90')
                self.write_f90(file_name)

    def write_f90(self, file_name=None, use='Cantera', routine_name='customkinetics',
                  mixture_name=None, fuel='FUEL', author='Unknown Author',
                  semi_implicit=False, rrate=True,
                  exponential=False, implicit_exponential=False, classical=True, dt_exp_user=1e-10,
                  supplement_dict=None, maintain_given_order=True, fit_order=None,
                  start_temperature=None, end_temperature=None, number_of_temperatures=None, pressure=None):
        """Writes the mechanism in f90 format

        Parameters
        ----------
        file_name :
            name of the output mechanism (Default value = None)
        use :
             (Default value = 'Cantera')
        routine_name :
             (Default value = 'customkinetics')
        mixture_name :
             (Default value = None)
        fuel :
             (Default value = 'FUEL')
        author :
             (Default value = 'Unknown Author')
        semi_implicit :
             (Default value = False)
        rrate :
             (Default value = True)
        exponential :
             (Default value = False)
        implicit_exponential :
             (Default value = False)
        classical :
             (Default value = True)
        dt_exp_user :
             (Default value = 1e-10)
        supplement_dict :
             (Default value = None)
        maintain_given_order :
             (Default value = True)
        fit_order :
             (Default value = None)
        start_temperature :
             (Default value = None)
        end_temperature :
             (Default value = None)
        number_of_temperatures :
             (Default value = None)
        pressure :
             (Default value = None)

        """

        if use == 'AVBP':
            if not mixture_name:

                warning_text = '\nWARNING : You should specify '
                plural_text = 0
                if fuel == 'FUEL':
                    warning_text += "a fuel name with the keyword fuel_name='FUEL' "
                    plural_text += 1
                if author == 'Unknown Author':
                    if 'ARCANE_AUTHOR' in os.environ:
                        author = os.environ.get('ARCANE_AUTHOR')
                        if plural_text == 1:
                            warning_text += '\nand '
                        warning_text += "an author name with the keyword author='Wade Wilson'"
                        plural_text += 1

                if plural_text > 0:
                    logger.warning(warning_text)

                # Creating initials from full_name
                initials_list = author.split(' ')
                if not len(initials_list) == 1:
                    initials = ''
                    for part in initials_list:
                        initials += part[0]
                else:
                    initials = author

                file_name = fuel + '_' + str(self.ns) \
                            + '_' + str(self.nr + self.nr_reverse) \
                            + '_' + str(self.nqss) + '_' + initials

                routine_name = file_name

            else:
                routine_name = mixture_name

        if not file_name:
            file_name = self.name + '.f90'

        if not fit_order:
            fit_order = self.fit_order

        if not start_temperature:
            start_temperature = self.start_temperature

        if not end_temperature:
            end_temperature = self.end_temperature

        if not number_of_temperatures:
            number_of_temperatures = self.number_of_temperatures

        if not pressure:
            pressure = self.fit_pressure

        custom.print_fortran(self, file_name, self.species_qss_names,
                             use=use, routine_name=routine_name,
                             semi_implicit=semi_implicit, rrate=rrate,
                             exponential=exponential, implicit_exponential=implicit_exponential,
                             classical=classical, dt_exp_user=dt_exp_user,
                             supplement_dict=supplement_dict, maintain_given_order=maintain_given_order,
                             fit_order=fit_order,
                             start_temperature=start_temperature,
                             end_temperature=end_temperature,
                             number_of_temperatures=number_of_temperatures,
                             pressure=pressure)

        logger.info(f".f90 file written as {file_name}")

    def write_chemkin(self, file_name=None, transport=True, thermo=True):
        """Writes the mechanism in Chemkin format

        Parameters
        ----------
        file_name :
            name of the output mechanism (Default value = None)
        transport :
            if True, writes a .tran file containing the transport data (Default value = True)
        thermo :
            if True, writes a .thermo file containing the thermodynamic data (Default value = True)

        """

        if not file_name:
            file_name = self.name + '.inp'

        if self.nqss:
            ct2ck.write(self.skeletal.ctmech, file_name, transport=transport, thermo=thermo)
        else:
            ct2ck.write(self.ctmech, file_name, transport=transport, thermo=thermo)

    def computed_cases(self, display=False, solution=False, create_case=True):
        """Generates a dictionary containing information about already computed cases

        Parameters
        ----------
        display :
            outputs the results of the search (Default value = False)
        solution :
            if True, looks for solution objects (Default value = False)
        create_case :
            if True, the class :func:`~ARCANE.cases.Case` object will be created and returned (Default value = True)
            else the ids are returned

        Returns
        -------
        meta_cases_list : list
            list of computed class :func:`~ARCANE.cases.Case` objects

        """

        computed_cases_id = database.get_cases_ids(self)

        if create_case:
            meta_cases_list = [tools.parse_id(case_id) for case_id in computed_cases_id]
            import ARCANE.cases as cases
            cases_list_temp = []
            for case in meta_cases_list:
                cases_list_temp.append(cases.create_case(reactor=case.reactor_type,
                                                         mechanism=self,
                                                         fuel=case.fuel,
                                                         oxidizer=case.oxidizer,
                                                         pressure=case.pressure,
                                                         temperature=case.temperature,
                                                         fuel_temperature=case.fuel_temperature,
                                                         phi=case.phi)[0])
            meta_cases_list = cases_list_temp

        else:
            meta_cases_list = computed_cases_id

        if display:
            for case in meta_cases_list:
                if create_case:
                    logger.info(case.myid)
                else:
                    logger.info(case)

        return meta_cases_list

    def specific_cases(self, reactor_type=None, P=-1, T=-1, phi=-1, fuel={}, oxidizer={}, display=True,
                       create_case=True, solution=False, clean=False):
        """Get a list of cases corresponding to the specific initial state given

        Parameters
        ----------
        reactor_type :
            type of reactor wanted (Default value = None)
        P :
            initial pressure (Default value = -1)
        T :
            initial temperature (Default value = -1)
        phi :
            equivalence ratio (Default value = -1)
        fuel :
            fuel directory (Default value = {})
        oxidizer :
            oxidizer directory (Default value = {})
        display :
            outputs the results of the search (Default value = True)
        solution :
            if True, looks for solution objects (Default value = False)
        create_case :
             (Default value = True)
        clean :
             (Default value = False)

        Returns
        -------
        specific_cases_list : list
            list of corresponding cases

        """

        meta_cases_list = self.computed_cases(solution=solution)#, create_case=create_case)
        specific_cases_list = tools.get_case_by_state(meta_cases_list, reactor_type=reactor_type,
                                                      P=P, T=T, phi=phi, fuel=fuel, oxidizer=oxidizer)

        # Sorting the data
        specific_cases_list = sorted(specific_cases_list, key=lambda x: (x.reactor_type,
                                                                         x.oxidizer,
                                                                         x.fuel,
                                                                         x.pressure,
                                                                         x.temperature,
                                                                         x.phi))

        if create_case:
            import ARCANE.cases as cases
            specific_cases_list_temp = []
            for case in specific_cases_list:
                specific_cases_list_temp.append(cases.create_case(reactor=case.reactor_type,
                                                                  mechanism=self,
                                                                  fuel=case.fuel,
                                                                  oxidizer=case.oxidizer,
                                                                  pressure=case.pressure,
                                                                  temperature=case.temperature,
                                                                  phi=case.phi)[0])
            specific_cases_list = specific_cases_list_temp

        if display:
            header_string = "\nHere is a list of available cases "
            if reactor_type:
                header_string += f"{reactor_type} "
            if P > 0:
                header_string += f"P = {P} Pa "
            if T > 0:
                header_string += f"T = {T} K "
            if phi > 0:
                header_string += f"phi = {phi} "
            if fuel:
                header_string += f"fuel: {fuel} "
            if oxidizer:
                header_string += f"oxidizer: {oxidizer} "

            logger.info(header_string)

            for case in specific_cases_list:
                info_string = ''
                if not reactor_type:
                    info_string += f"{case.reactor_type} "
                if P < 0:
                    info_string += f"P = {case.pressure} Pa "
                if T < 0:
                    info_string += f"T = {case.temperature} K "
                if phi < 0:
                    info_string += f"phi = {case.phi} "
                if not fuel:
                    info_string += f"fuel: {case.fuel} "
                if not oxidizer:
                    info_string += f"oxidizer: {case.oxidizer} "

                logger.info(info_string)

        if clean:
            logger.info("Beginning the cleaning of cases (Only available for 1D premixed flames)")
            for case in specific_cases_list:
                if case.reactor_type.startswith('C1D'):
                    logger.info(f"Cleaning {case.myid}")
                    case.clean_solutions(mechanism=self)

        return specific_cases_list


def cantera_reduced_solution(ctmech, species_to_discard=[], species_to_add=[],
                             reactions_to_discard=[], reactions_to_add=[],
                             diluent=[], auto_remove=False, reduced_cti_file=None,
                             transport=None):
    """Removes species or reactions from a Cantera Solution object

    Parameters
    ----------
    ctmech :
        Cantera Solution object
    species_to_discard :
        list of species to discard (Default value = [])
    species_to_add :
        list of species objects to add (Default value = [])
    reactions_to_discard :
        list of reactions indices to discard (Default value = [])
    reactions_to_add :
        list of reactions objects to add (Default value = [])
    diluent :
        species that even inactive will be kept (Default value = [])
    auto_remove :
        if True, removes species that are no longer used in a reaction (Default value = False)
    reduced_cti_file :
        name of the output .cti (Default value = None)
    transport :
        transport model (Default value = None)

    Returns
    -------
    phase : class `Cantera.Solution` object
        modified Cantera Solution object

    """

    # Species objects list
    species_objects = ctmech.species()

    # Species list
    species_names = ctmech.species_names

    # Reactions objects list
    reactions_objects = ctmech.reactions()

    # Species objects list sorted according to the species list
    species_objects_names = [spec.name for spec in species_objects]

    # Corrects inconsistency between the list of species and the list of species objects
    if species_objects_names != species_names:
        species_objects_sorted = [species_objects[species_objects_names.index(spec)] for spec in species_names]
    else:
        species_objects_sorted = species_objects

    # Species to include in new mechanism
    species_to_include = [s for s in species_objects_sorted if s.name not in species_to_discard]
    species_to_include_names = [s.name for s in species_to_include]

    species_in_a_reaction = set()

    # Initialization
    reactions_to_include = []

    # Reactions to include in new mechanism
    for reac_base in reactions_objects:
        add_reac = True

        if not all(reactant in species_to_include_names for reactant in reac_base.reactants):
            continue

        if not all(product in species_to_include_names for product in reac_base.products):
            continue

        if species_to_discard and reac_base.reaction_type in [2, 4]:
            new_efficiencies = reac_base.efficiencies
            species_in_a_reaction.update(new_efficiencies.keys())

            for spec in species_to_discard:
                if spec in reac_base.efficiencies:
                    new_efficiencies.pop(spec)
            reac_base.efficiencies = new_efficiencies

        for reac in reactions_to_discard:
            if tools.same_reaction(reac_base, reac):
                add_reac = False
                break

        if add_reac:
            # Species presents in a reaction
            species_in_a_reaction.update(reac_base.reactants)
            species_in_a_reaction.update(reac_base.products)

            reactions_to_include.append(reac_base)

    species_in_a_reaction.update(diluent)

    # Adding new reaction objects
    # TODO put a test on wether or not all the species are valid
    if reactions_to_add:
        reactions_to_include = reactions_to_add + reactions_to_include

    if species_to_add:
        species_to_include_names = [spec.name for spec in species_to_include]
        for spec in species_to_add:
            if spec.name not in species_to_include_names:
                species_to_include = [spec] + species_to_include

    # If a species is no more in a reaction, exclude it
    automatically_removed = list(set(species_to_include_names) - set(species_in_a_reaction))

    if automatically_removed and auto_remove:
        logger.debug('Species ' + ', '.join(automatically_removed)
                     + ' automatically removed since it was no longer present in a reaction')

        species_to_include = [s for s in species_to_include if s.name in species_in_a_reaction]

    # Generating the cti file (if qss_cti_file= None, a dummy one will be generated)
    if not reduced_cti_file:
        reduced_cti_file = "dummy.cti"

    phase = ct.Solution(thermo='IdealGas',
                        transport=transport,
                        kinetics='GasKinetics',
                        species=species_to_include,
                        reactions=reactions_to_include)

    ct2cti.write(phase, reduced_cti_file, transport=transport)

    if reduced_cti_file == "dummy.cti":
        subprocess.call('rm -rf ' + reduced_cti_file, shell=True)

    return phase
