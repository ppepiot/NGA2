"""
Functions handling the aesthetic aspects of the code
"""

import logging
import datetime
import time
import sys
import os

# Setting global variables
colours_on = True
cantera_solve_loglevel = 1
cantera_restore_loglevel = 0

level = 'info'
level_index = 2

colours = {}
colours['default'] = '\033[0m'

colours['bold'] = '\033[1m'
colours['underline'] = '\033[4m'
colours['overline'] = '\033[27m'

colours['black]'] = '\033[30m'
colours['grey'] = '\033[90m'

colours['red'] = '\033[91m'
colours['green'] = '\033[92m'
colours['yellow'] = '\033[93m'
colours['blue'] = '\033[94m'
colours['purple'] = '\033[95m'
colours['cyan'] = '\033[96m'

colours['red_2'] = '\033[31m'
colours['green_2'] = '\033[32m'
colours['yellow_2'] = '\033[33m'
colours['blue_2'] = '\033[34m'
colours['purple_2'] = '\033[35m'
colours['cyan_2'] = '\033[36m'

colours['highlight_black'] = '\033[40m'
colours['highlight_red'] = '\033[41m'
colours['highlight_green'] = '\033[42m'
colours['highlight_yellow'] = '\033[43m'
colours['highlight_blue'] = '\033[44m'
colours['highlight_purple'] = '\033[45m'
colours['highlight_cyan'] = '\033[46m'
    

class Logger(object):
    
    def __init__(self):
        """Overcharging existing class
        """ 

        self.colours_on = colours_on
        self.logger = self.set_log('logCompute')

        self.default_c = colours['default']
        self.debug_c = colours['grey']
        self.info_c = colours['default']
        self.goodnews_c = colours['green']
        self.warning_c = colours['yellow']
        self.error_c = colours['red']
        self.critical_c = colours['highlight_red']

        # Current logger level
        self.level = level
        self.level_index = level_index

    def terminator(self, term_string='\n'):
        """Sets the terminator character of the logger

        Parameters
        ----------
        term_string :
            string to end the logger (Default value = `newline`)

        """

        self.logger.handlers[1].terminator = term_string

    def debug(self, text):
        """Writes text in debug mode

        Parameters
        ----------
        text :
            text to write

        """

        if not isinstance(text, str):
            text = str(text)

        if colours_on:
            colour_string = self.debug_c
            end_string = self.default_c
            text = colour_string + text + end_string

        self.logger.debug(text)

    def info(self, text):
        """Writes text in information mode

        Parameters
        ----------
        text :
            text to write

        """

        if not isinstance(text, str):
            text = str(text)

        if colours_on:
            colour_string = self.info_c
            end_string = self.default_c
            text = colour_string + text + end_string

        self.logger.info(text)

    def goodnews(self, text):
        """Writes text for good news

        Parameters
        ----------
        text :
            text to write

        """

        if not isinstance(text, str):
            text = str(text)

        if colours_on:
            colour_string = self.goodnews_c
            end_string = self.default_c
            text = colour_string + text + end_string

        self.logger.info(text)

    def warning(self, text):
        """Writes text for warning

        Parameters
        ----------
        text :
            text to write

        """

        if not isinstance(text, str):
            text = str(text)

        if colours_on:
            colour_string = self.warning_c
            end_string = self.default_c
            text = colour_string + text + end_string

        self.logger.warning(text)

    def error(self, text):
        """Writes text for error

        Parameters
        ----------
        text :
            text to write

        """

        if not isinstance(text, str):
            text = str(text)

        if colours_on:
            colour_string = self.error_c
            end_string = self.default_c
            text = colour_string + text + end_string

        self.logger.error(text)

    def critical(self, text):
        """Writes text for critical message

        Parameters
        ----------
        text :
            text to write

        """

        if not isinstance(text, str):
            text = str(text)

        if colours_on:
            colour_string = self.critical_c
            end_string = self.default_c
            text = colour_string + text + end_string

        self.logger.critical(text)

    def set_log(self, log_name='logCompute', level=level, i=level_index):
        """Gets the correct logger

        Parameters
        ----------
        log_name :
            string for logger name (Default value = 'logCompute')
        level :
            name of the verbose level (debug, info, warning, error, critical) (Default value = 'info')
        i :
            corresponding number for the verbose level (1 to 5) (Default value = 2)

        Returns
        -------
        logger : class `ARCANE.logger` object
            class `ARCANE.logger` object to be handled

        """

        logger = logging.getLogger(log_name)

        if level == 'debug' or i == 1:
            logger.setLevel(logging.DEBUG)
            self.level = 'debug'
            self.level_index = 1

        if level == 'info' or i == 2:
            logger.setLevel(logging.INFO)
            self.level = 'info'
            self.level_index = 2

        if level == 'warning' or i == 3:
            logger.setLevel(logging.WARNING)
            print('set')
            self.level = 'warning'
            self.level_index = 3

        if level == 'error' or i == 4:
            logger.setLevel(logging.ERROR)
            self.level = 'error'
            self.level_index = 4

        if level == 'critical' or i == 5:
            logger.setLevel(logging.CRITICAL)
            self.level = 'critical'
            self.level_index = 5

        if level not in ['debug', 'info', 'warning', 'error', 'critical'] \
                and i not in [1, 2, 3, 4, 5]:
            logger.critical('Debug level will be set to info by default')
            logger.setLevel(logging.INFO)

        self.logger = logger
    
        return logger

    def create_log(self, name='File'):
        """Handles the verbose levels
        
        Parameters
        ----------
        name :
             name of the log file (Default value = 'File')

        """
        # Instantiate the logger object
        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)

        # Create the log file that can be read after the computation
        file_logger = logging.FileHandler(name + '.log', mode='w')
        logger.addHandler(file_logger)

        # Create the message occurring in the console
        console = logging.StreamHandler()
        logging.getLogger(name).addHandler(console)


internal_logger = Logger()
internal_logger.set_log('logCompute')


def display(status):
    """Display or not on terminal

    Parameters
    ----------
    status :
        on or off

    """

    if status in [True, "on", "yes"]:
        sys.stdout = sys.__stdout__
    elif status in [False, "off", "no"]:
        sys.stdout = open(os.devnull, 'w')
    else:
        internal_logger.warning('WARNING ! Argument not valid, code will carry on ...')


def print_colour(*object, colour='default'):
    """Prints the text in the desired colour

    Parameters
    ----------
    *object :
        object you want to print
    colour :
        colour (Default value = 'default')

    """

    if colour in colours and colours_on:
        print(colours[colour], end='')
        print(*object, end='')
        print(colours['default'])
    else:
        print(*object)


def ARCANE_logo(size='big'):
    """Printing the ARCANE logo
    
    Parameters
    ----------
    size :
         String to change the size of the ARCANE logo (Default value = 'big')

    """

    logo_big = """\
========================================================================

                                  ^                                  
                                 / \                                 
                                /   \                                
                               /     \                               
                              /      |\                            
                             /       | \                             
                            /        |  \                           
                           /        /    \                         
                          /        /      \                         
                         /        /        \                       
                        /      __/          \                       
                       /      /              \                    
                      /      /                \                    
                     / \    |                / \                     
                    /   \   |               /   \                 
                   /     \  |              /     \              
                  /    @  \  \            /  @    \           
                 /      \  \  \          /  /      \        
                /        \  \__\________/  /        \ 
                \     @---@       |   |   @---@     /
                 \    |    \  @-, | ,-@  /    |    /
                  \   @     \/   \|/   \/     @   /
                   \    \   @ ----@---- @   /    /
                    \    \   \    |    /   /    /
                     \    \   \   |   /   /    /
                      \    @    --@--    @    /
                       \  /  ____/ \____  \  /
                        \/  /           \  \/ 
                         -------------------

                     _____  _____            _    _ _____    
              /\    |     \/  ___|    /\    | \  | | ____|               
             /  \   | | > |  /       /  \   |  \ | | |___                  
            / /\ \  |  _  / |       / /\ \  | |\ \ |  ___|                   
           / ____ \ | | \ \  \___  / ____ \ | | \  | |___                 
          /_/    \_\|_|  \_\_____|/_/    \_\|_|  \_|_____|                       

              Property of CERFACS and Cornell University

========================================================================
"""

    logo_small = """\
========================================================================
                     _____  _____            _    _ _____
              /\    |     \/  ___|    /\    | \  | | ____|
             /  \   | | > |  /       /  \   |  \ | | |___
            / /\ \  |  _  / |       / /\ \  | |\ \ |  ___|
           / ____ \ | | \ \  \___  / ____ \ | | \  | |___
          /_/    \_\|_|  \_\_____|/_/    \_\|_|  \_|_____|

              Property of CERFACS and Cornell University

========================================================================
"""

    if size == 'small':
        logo = logo_small
    elif size == 'no':
        logo = ""
    else:
        logo = logo_big

    return logo


def banner(colour='auto'):
    """Returns a banner for terminal display

    Parameters
    ----------
    colour :
        colouration of the banner (Default value = 'auto')

    Returns
    -------
    banner : str
        brilliant banner art in a string

    """

    now = datetime.datetime.now()

    if (now.month == 12 and 20 < now.day < 31) \
            or colour == 'holidays':
        banner = \
            """
                .--._.--.--.__.--.--.__.--.--.__.--.
               _(_      _Y_      _Y_      _Y_      _)_
              [___]    [___]    [___]    [___]    [___]
              /:' \    /:' \    /:' \    /:' \    /:' \ 
             |::   |  |::   |  |::   |  |::   |  |::   |
             \::.  /  \::.  /  \::.  /  \::.  /  \::.  /
              \::./    \::./    \::./    \::./    \::./
               '='      '='      '='      '='      '='
            
                        HAPPY HOLIDAYS !!!
            
            """
    elif (now.month == 1 and now.day < 15) \
            or colour == 'new year':
        banner = \
            """
             O  .    O  .    O  .    O  .    O  .            
              .  O    .  O    .  O    .  O    .  O        
             O_._o   O_._o   O_._o   O_._o   O_._o           
            \~~~~~/ \~~~~~/ \~~~~~/ \~~~~~/ \~~~~~/     
             \   /   \   /   \   /   \   /   \   /     
              \ /     \ /     \ /     \ /     \ /         
               I       I       I       I       I        
               I       I       I       I       I        
             __I__   __I__   __I__   __I__   __I__      
            
                        HAPPY NEW YEAR !!!   
            
            """

    elif (now.month == 10 and now.day == 31) \
            or colour == 'halloween':
        banner = \
            """
               /\               /\ 
              /  \_   (\_/)   _/  \ 
             /     '--('.')--'     \ 
             |  _     /   \     _  |
              \/  \_|`\___/`|_/  \/
                     \(/ \)/       
                      "   "
            
            """
    elif now.month in [7, 8] or colour == 'summer':
        banner = \
            """
                     \ | /
                    --( )--
                     / | \       _\/_
                                 //o\   -\/_
               _____ _ __ __ ____ _ | __/o\_
             =-=-_-__=_-= _=_=-=_,-'|"'""-|-,_
              =- _=-=- -_=-=_,-"          |
                =- =- -=.--"
            
            """
    else:
        banner = ''

    return banner


def print_greetings_and_farewell(greeting=False, farewell=False):
    """Print greetings and farewell messages dependent of the time

    Parameters
    ----------
    greeting :
        if True prints the message (Default value = False)
    farewell :
        if True prints the message (Default value = False)

    Returns
    -------
    greeting_message : str
        greeting message at the beginning of the code
    farewell_message : str
        farewell message at the end of the code
        

    """

    # Error message

    now = datetime.datetime.now()

    greeting2_message = None
    if 7 <= now.hour < 12:
        greeting_message = 'Good morning,\n'
        farewell_message = 'Have a nice day.'
    elif 12 <= now.hour < 18:
        greeting_message = 'Good afternoon,\n'
        farewell_message = 'Have a nice day.'
    elif 18 <= now.hour < 23:
        greeting_message = 'Good evening,\n'
        farewell_message = 'Have a nice evening.'
    elif now.hour >= 23 or now.hour < 4:
        greeting_message = "Good ev... have you seen the time ?!\n" \
                           "Why are you working so late ?! Go to sleep, you'll do that tomorrow !"
        greeting2_message = "\nOkay you're still there ..."
        farewell_message = 'Go sleep now !!! Nighty night !'
    else:
        greeting_message = "Good morning,\nWhy are you awake so early ?! Go back to bed, dammit !"
        greeting2_message = "\nOkay you're still there ..."
        farewell_message = 'Have a nice (long) day (go take a coffee, you will need it) !'

    if greeting:
        print(greeting_message)
        if greeting2_message:
            time.sleep(5)
            greeting_message += greeting2_message
            print(greeting2_message)

    if farewell:
        print(farewell_message)

    return greeting_message, farewell_message


def head_tail_charac(text='', size='medium', style='#'):
    """Prints the headers for an output message

    Parameters
    ----------
    text :
        the text you want in your header (Default value = '')
    size :
        size of the header and tail (Default value = 'medium')
    style :
        character filling the header and tail (Default value = '#')

    Returns
    -------
    head_message : str
        head message string
    tail_message : str
        tail message string

    """

    half_line = style * 11
    head_text = half_line + ' ' + text + ' ' + half_line + '\n'
    line = style * (len(head_text) - 1) + '\n'

    if size == 'small':
        head_charac = head_text

        tail_charac = '\n' + line

    elif size == 'big':
        head_charac = (line * 2) + head_text + (line * 2)

        tail_charac = '\n' + line * 5

    else:
        head_charac = line + head_text + line

        tail_charac = '\n' + line * 3

    return head_charac, tail_charac
