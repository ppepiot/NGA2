"""Initialisation of the module

Author: QC (2018/12/19)
"""
import sys

version = sys.version.split('.')[1]

if int(version) <= 6:
    print('Your python version (anterior to 3.7) does not support circular import')
    print('You must import ARCANE.manual_init at the beginning of you script')
else:

    import ARCANE.database as database
    import ARCANE.display as display

    import atexit
    import os

    logger = display.Logger()

    # logger for automatic.py, drgep.py, qss_species.py, error.py
    logger.set_log('logReduction')
    logger.create_log('logReduction')

    # logger for the rest
    logger.set_log('logCompute')
    logger.create_log('logCompute')

    arcanoying = os.environ.get('ARCANE_GREETS')

    arcane_logo = os.environ.get('ARCANE_LOGO')

    arcane_colours = os.environ.get('ARCANE_COLOURS')

    if arcane_colours in ['off', 'no', 'n', 'nope', 'nein', 'niet']:
        arcane_colours = False
    else:
        arcane_colours = True

    display.colours_on = arcane_colours


    def exit_handler():
        """

        Handles what will be done at python exit

        :return: None
        """

        print()
        display.print_greetings_and_farewell(farewell=True)

        # Removing log files if they are empty
        number_of_lines = sum(1 for line in open('logCompute.log'))
        if number_of_lines <= 3:
            os.remove('logCompute.log')

        if not os.path.getsize('logReduction.log'):
            os.remove('logReduction.log')

        # Removing automatically generated folders if they are empty
        if system == 'database':
            if os.path.isdir('database') and not os.listdir('database'):
                os.rmdir('database')

        else:
            if os.path.isdir('cases') and not os.listdir('cases'):
                os.rmdir('cases')

            if os.path.isdir('mechanisms') and not os.listdir('mechanisms'):
                os.rmdir('mechanisms')


    if not arcanoying or arcanoying != 'silence':

        banner = display.banner()

        logo = display.ARCANE_logo(arcane_logo)

        display.print_colour(logo, colour='red')
        display.print_colour(banner)
        display.print_greetings_and_farewell(greeting=True)
        display.print_colour("Have fun using ARCANE !")
        print()

        atexit.register(exit_handler)

    else:

        logo = display.ARCANE_logo('small')
        display.print_colour(logo)

    # Initiating default database if ARCANE_DATABASE is a environment variable
    system = 'database'
    file_format = 'ascii'

    if 'ARCANE_DATA_SYSTEM' in os.environ:
        if os.environ['ARCANE_DATA_SYSTEM']:
            system = os.environ.get('ARCANE_DATA_SYSTEM')
            logger.debug('ARCANE_DATA_SYSTEM environment variable has been found.')

    if 'ARCANE_DATA_FORMAT' in os.environ:
        if os.environ['ARCANE_DATA_FORMAT']:
            file_format = os.environ.get('ARCANE_DATA_FORMAT')
            logger.debug('ARCANE_DATA_FORMAT environment variable has been found.')

    if 'ARCANE_DATABASE' in os.environ:
        if os.environ['ARCANE_DATABASE']:
            logger.debug('ARCANE_DATABASE environment variable has been found.')
            database.init_database('$ARCANE_DATABASE', system=system, file_format=file_format)
        else:
            database.init_database(path_to_database='.', system=system, file_format=file_format)
    else:
        database.init_database(path_to_database='.', system=system, file_format=file_format)
