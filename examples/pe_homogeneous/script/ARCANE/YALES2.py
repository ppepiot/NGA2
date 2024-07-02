"""Routines for direct integration into YALES2

Author: QC (2019/10/17)
"""

import os
import subprocess

import cantera.ctml_writer

import ARCANE.ct2cti as ct2cti
import ARCANE.display as display
import ARCANE.database as database
import ARCANE.graphs as graphs

logger = display.Logger()
logger.set_log('logCompute')

graphs = graphs.Graph()


def dump_YALES2_inputs(mechanism, output_dir="YALES2_input", mixture_name=None, exponential=False,
                       author='Unknown Author', since='2020.04', details='', note='', fuel_name='FUEL'):
    """Writes .xml and fortran files

    Parameters
    ----------
    mechanism : :func:`~ARCANE.mechanisms.Mechanism` object
        Object of class :func:`~ARCANE.mechanisms.Mechanism` object
    output_dir : str
        directory where the files will be written (Default value = "YALES2_input")
    mixture_name : str
        name of the mixture (Default value = None)
    exponential : bool
        if True uses the exponential formulation for production rates (Default value = False)
    author : str
        name of the mechanism author (for name of the mixture) (Default value = 'Unknown Author')
    since : str
        current version of YALES2 (useful to fill the f90 header) (Default value = '2020.04')
    details : str
        details on the f90 use (useful to fill the f90 header) (Default value = '')
    note : str
        note on reduction (useful to fill the f90 header) (Default value = '')
    fuel_name : str
        name of the fuel for automatic naming (Default value = 'FUEL')


    Created: 18/12/15 [QC]

    """

    # Create output directory if no already there (do not overwrite)
    database.create_dir(output_dir, overwrite=False)

    if not mixture_name:

        warning_text = '\nWARNING : You should specify '
        plural_text = 0
        if fuel_name == 'FUEL':
            warning_text += "a fuel name with the keyword fuel_name='FUEL' "
            plural_text += 1
        if author == 'Unknown Author':
            if 'ARCANE_AUTHOR' in os.environ:
                author = os.environ.get('ARCANE_AUTHOR')
            if plural_text == 1:
                warning_text += '\nand '
            warning_text += "an author name with the keyword author='TS' (if your name is Tony Stark)"
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

        mixture_name = fuel_name + '_' + str(mechanism.ns) \
                       + '_' + str(mechanism.nr + mechanism.nr_reverse) \
                       + '_' + str(mechanism.nqss) + '_' + initials

    if mechanism.species_qss_names:
        if mechanism.skeletal:
            mechanism_origin = mechanism.skeletal
        else:
            logger.info("No skeletal mechanism specified, the f90 file will not be written")
            return
    else:
        mechanism_origin = mechanism

    # Writes the .xml file
    file_name = mechanism.path.split('/')[-1].split('.')[0]
    cti_file = file_name + '.cti'

    ct2cti.write(mechanism.ctmech, cti_file, dummy="reaction('O => O', [0.0, 0.0, 0.0])")

    cantera.ctml_writer.convert(cti_file)
    subprocess.call('rm -rf ' + cti_file, shell=True)
    if os.path.isfile(file_name + '.xml'):
        subprocess.call('mv ' + file_name + '.xml ' + output_dir + '/' + mixture_name + '.xml', shell=True)
        logger.info(output_dir + '/' + mixture_name + '.xml successfully created')
    else:
        logger.error('ctml_writer encountered a problem, you need do create the xml file by yourself')

    # Header filling information
    supplementary_info = dict()
    supplementary_info['details'] = details
    supplementary_info['authors'] = author
    supplementary_info['since'] = since
    supplementary_info['note'] = note

    # Writes f90 file
    mechanism_origin.write_f90(file_name=mixture_name, use='YALES2', routine_name=mixture_name,
                               exponential=exponential, supplement_dict=supplementary_info)

    subprocess.call('mv ' + mixture_name + '.f90 ' + output_dir + '/' + mixture_name + '.f90', shell=True)
    logger.info(output_dir + '/' + mixture_name + '.f90 successfully created')
