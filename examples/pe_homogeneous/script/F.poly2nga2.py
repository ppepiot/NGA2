

from ckpoly2nga2 import *
from custom_printing_nga import *
from ARCANE import mechanisms
import numpy as np
from ARCANE import tools as tools

class ExtendedParser(Parser):

    def convertPoly(self,liquidFile, gasFile, thermoFile, permissive=False, ctiName=None):
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

        # NOTE: The following code is to deal with the third-body reactions
        #       But in Faravelli group's kinetic file, the third-body reactions are not standard.
        #       So this part of code can not be used in other places
        # for i in self.reactions:
        #     # find out thid-body reactions
        #     a = [(i.reactants[j][0],i.reactants[j][1].label) for j in range(len(i.reactants))]
        #     b = [(i.products[j][0],i.products[j][1].label) for j in range(len(i.products))]
        #     ThirdBodies = []
        #     for j in a:
        #         if j in b:
        #             ThirdBodies.append(j[1])
        #     if len(ThirdBodies) > 0:
        #         efficiencies = {}
        #         for j in ThirdBodies:
        #             efficiencies[j] = 1.0
        #         tbkinetic = ThirdBody(arrheniusHigh=i.kinetics, 
        #                             parser=self,
        #                             efficiencies=efficiencies)
        #         i.kinetics = tbkinetic
        #         i.thirdBody = ThirdBodies

        for i in self.reactions:
            # check duplicate reactions
            if i.duplicate:
                really_duplicate = False
                for j in self.reactions:
                    if i.reactants == j.reactants and i.products == j.products and i!=j:
                        really_duplicate = True
                i.duplicate = really_duplicate

        if ctiName is None:
            ctiName = os.path.splitext(liquidFile)[0] + '.cti'
        
        # Write output file
        surface_names = self.writeCTI(name="pyrolysis", outName=ctiName)
        print('Wrote CTI mechanism file to {0!r}.'.format(ctiName))
        print('Mechanism contains {0} species and {1} reactions.'.format(len(self.speciesList), len(self.reactions)))
        return surface_names

class react():

    def __init__(self):
        self.reactants = {}
        self.products = {}
        self.orders = {}

def modify_species_name(name):
    if '-' in name:
        name = name.replace('-','X')
    if '(' in name:
        name = name.replace('(','G')
    if ')' in name:
        name = name.replace(')','G')
    return name

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
        ctiName = outName[:-4] + '.cti'
    else:
        outName = os.path.splitext(liquidFile)[0] + '.f90'
        ctiName = os.path.splitext(liquidFile)[0] + '.cti'

    parser.convertPoly(liquidFile, gasFile, thermoFile,permissive=True, ctiName=ctiName)

    print("\n\n\n\n\n\n")

    # for reac in parser.reactions:
    #     a = [(reac.reactants[i][0],reac.reactants[i][1].label) for i in range(len(reac.reactants))]
    #     b = [(reac.products[i][0],reac.products[i][1].label) for i in range(len(reac.products))]
    #     for i in a:
    #         if i in b:
    #             reac.thirdBody = True

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

    use = "NGA2"

    print_header(print_variables, use)

    # prepare infomation for declaration
    print_variables['precision'] = 'WP'

    constants = {'Rcst': 8.314,'ns': 0,'nr': 0,'nr_reverse': 0,'nFO': 0,
    'nFO_reverse': 0,'nTB': 0,'nTB_reverse': 0,'nPlog': 0,'nPlog_reverse': 0,'nqss': 0,'ne': 0}
    constants['ns']=len(parser.speciesList)
    constants['nr']=len(parser.reactions)
    constants['ne']=len(parser.elements)

    # PAY ATTENTION : as in Faravelli's group, the pyrolysis has no fall-off reaction. We ignore it here.
    
    mech_variables = {}
    mech_variables['elements'] = parser.elements
    mech_variables['mech'] = type('mech', (object,), {'molecular_weights': [0.0]*len(parser.speciesList)})
    mech_variables['species'] = [parser.speciesList[i].label for i in range(len(parser.speciesList))]

    for i in range(len(parser.speciesList)):
        #calculate molecular weight
        weights = {'H':1.00794,'C':12.011,'N':14.00674,'O':15.9994,'He':4.002602}
        mw = 0
        composition = parser.speciesList[i].composition
        for j in list(composition.keys()):
            if j not in weights.keys():
                print('Error: Element '+j+' is not in the list of weights')
            else:
                mw += weights[j]*composition[j]
        mech_variables['mech'].molecular_weights[i] = mw
        if 'P' in parser.speciesList[i].label:
            mech_variables['mech'].molecular_weights[i] = -1.0
    
    # Correct species names in f90 
    mech_variables['species'] = tools.convert_to_valid_fortran_var(mech_variables['species'])

    qss_variables = {
        'species_qss_names':[]
    }

    reactions_variables = {
        'reac_label': [],
        'reac_direction': [],
        'reac_names': [],
        'ThreeBodiesindex': [],
        'FallOffindex': [],
        'Plogindex': [],
        'Plog_pressures': [],
    }

    count = 0
    for i in parser.reactions:
        reac_label = count
        reac_direction = ''
        if i.reversible:
            reac_direction = '_r'
        else:
            reac_direction = '_f'
        reactions_variables['reac_label'].append(reac_label)
        reactions_variables['reac_direction'].append(reac_direction)
        reactions_variables['reac_names'].append(str(i))
        # if i.thirdBody:
        #     reactions_variables['ThreeBodiesindex'].append(count)
        count += 1

    print_declarations(print_variables, constants, mech_variables, qss_variables, reactions_variables, use, semi_implicit_bool=False)
    
    print_species_names(print_variables,use, mech_variables)
    
    print_reaction_expressions(print_variables, constants,reactions_variables)

    # Following 3 functions are not used
    # print_lindemann(print_variables)

    # print_third_body(print_variables, constants, mech_variables, reactions_variables)

    arrhenius_constants={}
    arrhenius_constants['forward_arrh_pdep'] = []
    arrhenius_constants['reverse_arrh_pdep'] = []

    reactions_variables['reac_type'] = [1]*len(parser.reactions)

    # print_pdep(print_variables, constants, arrhenius_constants, mech_variables, reactions_variables)

    arrhenius_constants['forward_arrh'] = np.zeros((3,len(parser.reactions)))
    arrhenius_constants['forward_arrh_0'] = []
    arrhenius_constants['forward_arrh_inf'] = []

    arrhenius_constants['reverse_arrh'] = []
    arrhenius_constants['reverse_arrh_0'] = []
    arrhenius_constants['reverse_arrh_inf'] = []

    for k,i in enumerate(parser.reactions):
        arrhenius_constants['forward_arrh'][0,k]=i.kinetics.A[0]
        arrhenius_constants['forward_arrh'][1,k]=i.kinetics.b
        arrhenius_constants['forward_arrh'][2,k]=i.kinetics.Ea[0]

    print_rate_coeff(print_variables, constants, arrhenius_constants, mech_variables, reactions_variables)

    mech_variables['nur'] = np.zeros(len(parser.reactions))
    mech_variables['nup'] = np.zeros(len(parser.reactions))
    mech = mech_variables['mech']
    mech.reaction = []

    for i in range(len(parser.reactions)):
        reac = parser.reactions[i]
        mech_variables['nur'][i] = len(reac.reactants)
        mech_variables['nup'][i] = len(reac.products)
        reaction = react()
        for j in range(len(reac.reactants)):
            reaction.reactants[reac.reactants[j][1].label] = reac.reactants[j][0]
        for j in range(len(reac.products)):
            reaction.products[reac.products[j][1].label] = reac.products[j][0]
        reaction.orders = reac.fwdOrders
        # convert all keys into valid fortran variables
        kys = list(reaction.orders.keys())
        for j in kys:
            if j == convert2fvar(j):
                reaction.orders[j] = float(reaction.orders[j])
            else:
                reaction.orders[convert2fvar(j)] = float(reaction.orders[j])
                del reaction.orders[j]
        kys = list(reaction.reactants.keys())
        for j in kys:
            if j == convert2fvar(j):
                pass
            else:
                reaction.reactants[convert2fvar(j)] = float(reaction.reactants[j])
                del reaction.reactants[j]
        kys = list(reaction.products.keys())
        for j in kys:
            if j == convert2fvar(j):
                pass
            else:
                reaction.products[convert2fvar(j)] = float(reaction.products[j])
                del reaction.products[j]
        mech.reaction.append(reaction)

    print_reaction_rates(print_variables, constants, mech_variables, qss_variables, reactions_variables)

    print_variables['EM_flag'] = True

    print_reaction_rate_rhs_Jac(print_variables,constants, mech_variables,reactions_variables)

    print_mass_to_concentration(print_variables,use)

    print_variables['f'].close()



if __name__ == '__main__':
    main(sys.argv[1:])
    print("Done")