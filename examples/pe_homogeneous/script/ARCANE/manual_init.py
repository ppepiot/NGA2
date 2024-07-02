import ARCANE.tools as tools
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
    """Handles what will be done at python exit

    """

    print()
    display.print_greetings_and_farewell(farewell=True)


if not arcanoying or arcanoying != 'silence':

    banner = display.banner()

    logo = display.ARCANE_logo(arcane_logo)

    display.print_colour(banner)
    display.print_colour(logo, colour='red')
    display.print_greetings_and_farewell(greeting=True)
    display.print_colour("Have fun using ARCANE !")
    print()

    atexit.register(exit_handler)

else:

    logo = display.ARCANE_logo('small')
    display.print_colour(logo)

# Initiating default database if ARCANE_DATABASE is a environment variable
system = 'database'

if 'ARCANE_DATA_SYSTEM' in os.environ:
    if os.environ['ARCANE_DATA_SYSTEM']:
        system = os.environ.get('ARCANE_DATA_SYSTEM')
        logger.debug('ARCANE_DATA_SYSTEM environment variable has been found.')

if 'ARCANE_DATABASE' in os.environ:
    if os.environ['ARCANE_DATABASE']:
        logger.debug('ARCANE_DATABASE environment variable has been found.')
        tools.init_database('$ARCANE_DATABASE', system=system)
    else:
        tools.init_database(path_to_database='.', system=system)
else:
    tools.init_database(path_to_database='.', system=system)