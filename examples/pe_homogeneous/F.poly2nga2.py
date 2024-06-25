

from ckpoly2nga2 import *
from custom_printing_nga import *

class ExtendedParser(Parser):

    def convertPoly(self,liquidFile, gasFile, thermoFile, permissive=False):
        if liquidFile:
            liquidFile = os.path.expanduser(liquidFile)
        if gasFile:
            gasFile = os.path.expanduser(gasFile)
        if thermoFile:
            thermoFile = os.path.expanduser(thermoFile)
        if permissive is not None:
            self.warning_as_error = not permissive
        logging.basicConfig(level=logging.INFO)

        if not os.path.exists(gasFile) or not os.path.isfile(liquidFile):
            raise InputParseError('Liquid/Gas file not found: ' + liquidFile + ' ' + gasFile)
        try:
            self.loadChemkinFile(gasFile)
        except Exception as e:
            logging.warning("\nERROR: Unable to parse '{0}' near line {1}:\n".format(
                                gasFile, self.line_number))
            raise
        try:
            self.loadChemkinFile(liquidFile)
        except Exception as e:
            logging.warning("\nERROR: Unable to parse '{0}' near line {1}:\n".format(
                                liquidFile, self.line_number))
            raise

        if not os.path.exists(thermoFile):
            raise IOError('Missing thermo file: {0!r}'.format(thermoFile))
        try:
            self.loadChemkinFile(thermoFile,skipUndeclaredSpecies=bool(liquidFile))
        except Exception:
            logging.warning("\nERROR: Unable to parse '{0}' near line {1}:\n".format(thermoFile, self.line_number))
            raise
        

def main(argv):

    print("\n\n\n\n\n\n")

    longOptions = ['liquid=', 'thermo=', 'transport=', 'gas=', 'id=',
                   'output=', 'permissive', 'help', 'debug']

    try:
        optlist, args = getopt.getopt(argv, 'dh', longOptions)
        options = dict()
        for o,a in optlist:
            options[o] = a

        if args:
            raise getopt.GetoptError('Unexpected command line option: ' +
                                     repr(' '.join(args)))

    except getopt.GetoptError as e:
        print('ck2cti.py: Error parsing arguments:')
        print(e)
        print('Run "ck2cti.py --help" to see usage help.')
        sys.exit(1)

    parser = ExtendedParser()

    if not options or '-h' in options or '--help' in options:
        print("""
This file is to transfer polymer pyrolysis kinetic file from chemkin format to NGA2(Fortran) format.
Usage:
    --liquid=*.liquid
    --gas=*.gas
    --thermo=*.tdc
    --output=*
              """)
        sys.exit(0)

    if '--liquid' in options:
        liquidFile = options['--liquid']
    else:
        liquidFile = None
    
    if '--gas' in options:
        gasFile = options['--gas']
    else:
        gasFile = None
    
    if '--thermo' in options:
        thermoFile = options['--thermo']
    else:
        thermoFile = None
    
    if (liquidFile is None) or (gasFile is None) or (thermoFile is None):
        print('Error: No input file specified.')
        sys.exit(1)

    if '--output' in options:
        outName = options['--output']
        if not outName.endswith('.f90'):
            outName += '.f90'
    else:
        outName = os.path.splitext(liquidFile)[0] + '.f90'

    permissive = '--permissive' in options

    parser.convertPoly(liquidFile, gasFile, thermoFile,permissive=permissive)

    print("\n\n\n\n\n\n")

    for reac in parser.reactions:
        a = [(reac.reactants[i][0],reac.reactants[i][1].label) for i in range(len(reac.reactants))]
        b = [(reac.products[i][0],reac.products[i][1].label) for i in range(len(reac.products))]
        for i in a:
            if i in b:
                reac.thirdBody = True

    # prepare infomation for header
    print_variables = {
        'f': open(outName, 'w'), 
        'f90_filename': outName,
        'routine_name': 'pyrolysis',
        'rrate': True,
        'exponential': False,
        'supp_info': {
            'details': 'None',
            'authors': 'None',
            'since': 'None',
            'note': 'This is a generated file.'
        }
    }
    use = 'Cantera'

    print_header(print_variables, use)

    # prepare infomation for declaration
    print_variables['precision'] = 'WP'

    constants = {'Rcst': 8.314,'ns': 0,'nr': 0,'nr_reverse': 0,'nFO': 0,
    'nFO_reverse': 0,'nTB': 0,'nTB_reverse': 0,'nPlog': 0,'nPlog_reverse': 0,'nqss': 0,'ne': 0}
    constants['ns']=len(parser.speciesList)
    constants['nr']=len(parser.reactions)
    constants['ne']=len(parser.elements)
    for i in parser.reactions:
        if i.reversible:
            constants['nr_reverse']+=1
        if i.thirdBody:
            constants['nTB']+=1
    
    # PAY ATTENTION : as in Faravelli's group, the pyrolysis has no fall-off reaction. We ignore it here.
    
    mech_variables = {}
    mech_variables['elements'] = parser.elements
    mech_variables['mech'] = type('mech', (object,), {'molecular_weights': []})
    mech_variables['species'] = [parser.speciesList[i].label for i in range(len(parser.speciesList))]
    mech_variables['mech'].molecular_weights = [0 for i in range(len(parser.speciesList))]

    qss_variables = {}

    # 示例的 reactions_variables 字典
    reactions_variables = {
        'reac_label': [1, 2, 3, 4, 5],
        'reac_direction': ['+', '-', '+', '-', '+'],
        'ThreeBodiesindex': [0, 2, 4],
        'FallOffindex': [1, 3],
        'Plogindex': [0, 2],
        'Plog_pressures': [2, 2, 2, 2, 2]
    }

    print_declarations(print_variables, constants, mech_variables, qss_variables, reactions_variables, use, semi_implicit_bool=False)

    print_variables['f'].close()



if __name__ == '__main__':
    main(sys.argv[1:])
    print("Done")