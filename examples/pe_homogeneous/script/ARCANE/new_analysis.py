import numpy as np
import os
#Graphs
import networkx as nx
import graphviz as gv
#ARCANE
import ARCANE.fileio as fileio
import ARCANE.sampling as sampling
import ARCANE.postproc as postproc
import ARCANE.display as display
#Logging
logger = display.Logger()
logger.set_log('logAnalysis')
logger.create_log('logAnalysis')

# Dumping
global analysis_dir
analysis_dir = 'analysis_results'
fileio.create_dir(analysis_dir, overwrite=False)

class Matrices:
    def __init__(self, caselist, mechanism, analysis, method, atom='C', resize_method=None, resize_parameter=None, n_samples=10):
        self.caselist = caselist
        self.mechanism = mechanism
        self.analysis = analysis
        self.method = method
        self.atom = atom
        self.resize_method = resize_method
        self.resize_parameter = resize_parameter
        self.n_samples = n_samples

        self.matrix, self.matrix_integrated, self.data_dict = self.caselist_to_paths()

    def caselist_to_paths(self):
        # Create sample database on cases and set-up solution file
        logger.info("=> Creating sampling database")

        # Build the database of databases to compare 2 cases
        # Maybe add a loop to calculate several mechanisms

        # data_dict[nvar] = [T1, T2, ..., Tnsamples]
        # MATRIX[case] =
        # DATADICT[case]['Y']
        # data_dict['T'] = [......]
        # MATRIX[ns, ns, sample, case], DATADICT[case] = create_matrix(case_list, resize_parameter, etc ...)
        sdbss = []
        matrix = {}
        matrix_integrated = {}
        data_dict = {}

        for nb_case, case in enumerate(self.caselist):
            resize_parameter_case = [self.resize_parameter[nb_case] if self.resize_parameter is list else self.resize_parameter]
            data_dict_case = self.resize_sdb(case, resize_parameter_case)
            data_to_sample_list, samples_database, remake = sampling.preprocess_samples_database(['T', 'P', 'Y', 'HR'],
                                                                                                self.mechanism, False)
            data_dict_final = sampling.extract_samples(data_dict_case, data_to_sample_list, integrated=False, n_samples=self.n_samples)
            sdbcase = sampling.postprocess_samples_database(self.mechanism, data_dict_final, remake)
            sdbss.append(sdbcase)
            data_dict[case] = data_dict_final

        # Loop on sdbss to analyse each case
        for index_sdb, sdb in enumerate(sdbss):

            # Compute linking matrix
            logger.info("=> Computing link matrix")

            if self.method == 'PFA':
                matrix_spec, matrix_reac = self.pfa_path(sdb)
            elif self.method == 'ATOM':
                matrix_spec, matrix_reac = self.atom_path(sdb)
            else:
                matrix_spec = []
                matrix_reac = []
                logger.error('Your method is not implemented for the moment !')
            
            if self.analysis == 'Pol':
                matrix_spec = np.transpose(matrix_spec, (1, 0, 2))

            if self.analysis in ['StoS', 'Glob', 'Pol']:
                matrix[self.caselist[index_sdb]] = matrix_spec
            else:
                matrix[self.caselist[index_sdb]] = matrix_reac

            # Compute integral matrix
            matrix_integrated[self.caselist[index_sdb]] = np.trapz(matrix[self.caselist[index_sdb]], sdb['grid'], axis=2)

        return matrix, matrix_integrated, data_dict

    def resize_sdb(self, case, resize_parameter_case):
        """ Function to resize the sample database.

        :param case: Case.
        :param resize_method: Resize method (string)
        :param resize_parameter_case: Resize paramater for the given case (array)
        :return: case
        """

        data_dict = case.data_dict(case.mechanism)
        if self.resize_method == 'Stencil':
            i_start = np.argmin(abs(data_dict['grid'] - resize_parameter_case[0]))
            i_end = np.argmin(abs(data_dict['grid'] - resize_parameter_case[1]))
            for key, value in data_dict.items():
                if key != 'mechanism':
                    data_dict[key] = value[i_start:i_end]
        elif self.resize_method == 'Point':
            i_point = np.argmin(abs(data_dict['grid'] - resize_parameter_case))
            for key, value in data_dict.items():
                if key != 'mechanism':
                    data_dict[key] = value[i_point]
        elif self.resize_method == 'HRR':
            # Definition in Rochette PhD multiplied by two
            x = data_dict['Grid']
            HRR = data_dict['HR']
            dHRR = np.gradient(HRR, x)
            d2HRR = np.gradient(dHRR, x)
            i_max = np.argmax(HRR)
            x_start = 2* (-dHRR[i_max] + np.sqrt(dHRR[i_max]**2 - 2*d2HRR[i_max]*HRR[i_max]))/d2HRR[i_max]
            x_end = 2* (-dHRR[i_max] - np.sqrt(dHRR[i_max]**2 - 2*d2HRR[i_max]*HRR[i_max]))/d2HRR[i_max]
            i_start = np.argmin(abs(x - (x[i_max]+x_start)))
            i_end = np.argmin(abs(x - (x[i_max]+x_end)))
            for key, value in data_dict.items():
                if key != 'mechanism':
                    data_dict[key] = value[i_start:i_end]

        return data_dict

    def pfa_path (self, sdb):
        """

        :param mechanism:
        :param data_dict:
        :return:
        """

        # Parameters
        n_data = len(sdb['grid'])
        n_spec = self.mechanism.ns
        n_reac = self.mechanism.nr
        net = self.mechanism.network
        nup = net.nup
        nur = net.nur
        nu = nup - nur

        sdb = postproc.extend_data_dict(sdb)

        DIC_spec = np.zeros([n_spec, n_spec, n_data], 'd')

        # Evaluate DIC(i,j)
        for i in range(n_spec):
            # reactions containing species i
            booli = net.indr[net.inds == i]
            indj = net.indsj[net.indsi == i]

            for j in indj:
                boolj = net.indr[net.inds == j]  # reactions containing species j
                indk = np.intersect1d(booli, boolj)  # reactions containing species i and j

                # Compute the DIC
                for k in indk:
                    if np.sign(nu[i, k]) != np.sign(nu[j, k]): # here is the modification
                        DIC_spec[i, j, :] += (nup[i, k] - nur[i, k]) * sdb['net rates of progress'][k, :]

                # Normalize
                DIC_spec[i, j, :] = DIC_spec[i, j, :] / np.maximum(np.maximum(sdb['creation rates'][i, :],
                                                                sdb['destruction rates'][i, :]), 1e-60)

        DIC_reac = np.zeros([n_spec, n_reac, n_data], 'd')

        # Evaluate DIC(i,j)
        for i in range(n_spec):
            # reactions containing species i
            indk = net.indr[net.inds == i]  # reactions containing species i

            # Compute the DIC
            for k in indk:
                DIC_reac[i, k, :] = (nup[i, k] - nur[i, k]) * sdb['net rates of progress'][k, :]

                # Normalize
                DIC_reac[i, k, :] = DIC_reac[i, k, :] / np.maximum(np.maximum(sdb['creation rates'][i, :],
                                                            sdb['destruction rates'][i, :]), 1e-60)

        return DIC_spec, DIC_reac

    def atom_path (self, sdb):
        """

        :param mechanism:
        :param data_dict:
        :return:
        """

        # Parameters
        n_data = len(sdb['grid'])
        n_spec = self.mechanism.ns
        ctmech = self.mechanism.ctmech
        n_reac = self.mechanism.nr
        net = self.mechanism.network
        nup = net.nup
        nur = net.nur
        nu = nup - nur

        sdb = postproc.extend_data_dict(sdb)
        omega = sdb['net rates of progress']

        # Workspace
        # DIC_spec = np.zeros([n_spec + 1, n_spec, n_data], 'd')
        DIC_spec = np.zeros([n_spec, n_spec, n_data], 'd')
        DIC_reac = np.zeros([n_spec + 1, n_reac, n_data], 'd')

        # Evaluate DIC(i,j)
        for i in range(n_spec):
            # reactions containing species i
            booli = net.indr[net.inds == i]
            indj = net.indsj[net.indsi == i]

            for j in indj:
                boolj = net.indr[net.inds == j]  # reactions containing species j
                indk = np.intersect1d(booli, boolj)  # reactions containing species i and j

                # Compute the reaction path A_ijk = (nA_i*nA_j*r_k)/NA_k
                for k in indk:
                    N = 0
                    for r in ctmech.reaction(k).reactants:
                        index = ctmech.species_index(r)
                        N += ctmech.n_atoms(r, self.atom) * nur[index, k]

                    signomega = np.sign(omega[k, :])
                    if N != 0:
                        forwardDIC = ctmech.n_atoms(ctmech.species_name(i), self.atom) * nur[i, k] * \
                                            ctmech.n_atoms(ctmech.species_name(j), self.atom) * nup[j, k] * omega[k, :] / N
                        backwardDIC = ctmech.n_atoms(ctmech.species_name(i), self.atom) * nur[j, k] * \
                                            ctmech.n_atoms(ctmech.species_name(j), self.atom) * nup[i, k] * omega[k, :] / N
                        DIC_spec[i, j, :] += forwardDIC * np.maximum(signomega, 0) \
                                                + backwardDIC * np.minimum(signomega, 0)

        return DIC_spec, DIC_reac
    
class Graph:
    def __init__(self, matrices, target_species, threshold=1e-1, size='FullHD', spec2plot=None):
        self.matrices = matrices
        self.target_species = target_species
        self.target_species_idx = [self.matrices.mechanism.ctmech.species_index(spec) for spec in self.target_species]
        self.threshold = threshold
        self.size_vector, self.ratio = self.get_size(size)
        self.spec_label = self.get_spec_label()
        self.reac_label = self.get_reac_label()
        self.spec2plot, self.spec2plot_idx = self.get_spec2plot(spec2plot)

        self.digraphs, self.digraphs_integrated = self.matrices_to_graphs()

    def get_spec2plot(self, spec2plot):
        print(self.matrices.mechanism.ctmech.species_names)
        if spec2plot is None:
            spec2plot = self.matrices.mechanism.ctmech.species_names
        elif isinstance(spec2plot, str):
            if '<=' in spec2plot:
                spec2plot = [s for s in self.matrices.mechanism.ctmech.species_names if self.matrices.mechanism.ctmech.n_atoms(s,spec2plot[0]) <= int(spec2plot[3:])]
            elif'==' in spec2plot:
                spec2plot = [s for s in self.matrices.mechanism.ctmech.species_names if self.matrices.mechanism.ctmech.n_atoms(s,spec2plot[0]) == int(spec2plot[3:])]
            elif '>=' in spec2plot:
                spec2plot = [s for s in self.matrices.mechanism.ctmech.species_names if self.matrices.mechanism.ctmech.n_atoms(s,spec2plot[0]) >= int(spec2plot[3:])]
            elif '<' in spec2plot:
                spec2plot = [s for s in self.matrices.mechanism.ctmech.species_names if self.matrices.mechanism.ctmech.n_atoms(s,spec2plot[0]) < int(spec2plot[2:])]
            elif '>' in spec2plot:
                spec2plot = [s for s in self.matrices.mechanism.ctmech.species_names if self.matrices.mechanism.ctmech.n_atoms(s,spec2plot[0]) > int(spec2plot[2:])]
        elif not isinstance(spec2plot, list):
            raise TypeError('spec2plot argument must be either list of species or condition on atom number "C<3"')
        
        spec2plot_idx = [self.matrices.mechanism.ctmech.species_index(spec) for spec in spec2plot]

        return spec2plot, spec2plot_idx

    def get_spec_label(self):
        spec_label = {}
        for indx in range(0, self.matrices.mechanism.ns):
            spec_label[indx] = str(self.matrices.mechanism.ctmech.species_name(indx))
        return spec_label

    def get_reac_label(self):
        reac_label = {}
        for indx in range(0, self.matrices.mechanism.nr):
            reac_label[indx] = str(self.matrices.mechanism.ctmech.reaction(indx))
        return reac_label

    def get_size(self, size):
        self.importance = 0.4
        # Check size of the plot
        if size == 'QHD':
            size_vector = [2560, 1440]
        elif size == 'FullHD':
            size_vector = [1920, 1080]
        elif size == 'HD':
            size_vector = [1280, 720]
        elif size == 'qHD':
            size_vector = [960, 540]
        elif size == 'YouTube144p':
            size_vector = [256, 144]
        elif isinstance(size, list):
            size_vector = [size[0], size[1]]
        else:
            logger.error('ERROR ! You should specify either the size of the vector or the keywords '
                        'QHD, FullHD, HD, qHD, Youtube144p.')
            sys.exit()
        ratio = size_vector[1] / ((1 - self.importance) * size_vector[0])
        return size_vector, ratio

    def matrices_to_graphs(self):
        digraphs = {}
        digraphs_integrated = {}
        for nb_case, case in enumerate(self.matrices.caselist):
            digraphs[case] = {}
            for i_sample, grid in enumerate(self.matrices.data_dict[case]['grid'][:]):
                matrix2plot = self.matrices.matrix[case][:, :, i_sample]
                # Remove unwanted species (to be improved)
                for i in range(len(matrix2plot)):
                    if i not in self.spec2plot_idx:
                        matrix2plot[i, :] = 0
                        matrix2plot[:, i] = 0
                digraphs[case][i_sample] = self.matrix_to_graph(matrix2plot, percentage=False)

            matrix2plot = self.matrices.matrix_integrated[case][:, :]
            # Remove unwanted species (to be improved)
            for i in range(len(matrix2plot)):
                if i not in self.spec2plot_idx:
                    matrix2plot[i, :] = 0
                    matrix2plot[:, i] = 0
            digraphs_integrated[case] = self.matrix_to_graph(matrix2plot, percentage=True)

        
        self.digraphs, self.digraphs_integrated = digraphs, digraphs_integrated
        return digraphs, digraphs_integrated

    def matrix_to_graph(self, matrix, percentage):
        logger.info("=> Drawing graph for links")

        if percentage:
            for i in range(len(matrix[:, 0])):
                sum_fluxes = np.sum(matrix[i, :])
                if sum_fluxes > 0.0:
                    for j in range(len(matrix[0, :])):
                        matrix[i, j] = matrix[i, j] / sum_fluxes * 100
            maximum = 100
        else:
            maximum = np.amax(abs(matrix))
        # Build digraph
        if self.matrices.analysis == 'StoRtoS':
            digraph = self.reactions_graph(matrix, percentage, maximum)
        elif self.matrices.analysis == 'RtoS':
            digraph = self.reactions_graph(matrix, percentage, maximum)
        elif self.matrices.analysis == 'StoS':
            digraph = self.species_graph(matrix, percentage, maximum)
        elif self.matrices.analysis in ['Glob', 'Pol']:
            digraph = self.species_graph(matrix, percentage, maximum)
        return digraph

    def species_graph(self, matrix, percentage, maximum):
        '''
        Return graph for StoS, Glob or Pol analysis
        :param matrix: Matrix with link between species
        :param maxim: Maximum of the DIC over all samples
        :param threshold: Minimum link weight to keep
        
        :return digraph: Digraph object
        '''
        # Filter DIC
        filtered_matrix = self.filter_species_matrix(matrix, maximum)
        if self.matrices.analysis == 'Pol':
            filtered_matrix = np.transpose(filtered_matrix)
        # Create graph
        digraph = nx.from_numpy_array(filtered_matrix, create_using=nx.DiGraph)
        #Remove nodes with no links
        isolated_nodes = list(nx.isolates(digraph))
        isolated_nodes = [node for node in isolated_nodes if node not in self.target_species_idx]
        digraph.remove_nodes_from(isolated_nodes)
        # Label nodes as species
        digraph = nx.relabel_nodes(digraph, self.spec_label)
        # Prepare graph for rendering
        self.tune_graph(digraph, percentage, maximum)

        return digraph

    def reactions_graph(self, matrix, percentage, maximum):
        filtered_matrix = self.filter_reactions_matrix(matrix, maximum)
        digraph = nx.DiGraph()
        for spec_idx in self.target_species_idx:
            digraph.add_node('spec%d' % spec_idx, label=self.spec_label[spec_idx], bipartite='species')
            for reac_idx in self.reac_label:
                if filtered_matrix[spec_idx, reac_idx] != 0:
                    digraph.add_node('reac%d' % reac_idx, label=self.reac_label[reac_idx], bipartite='reactions')
                    if filtered_matrix[spec_idx, reac_idx] > 0:
                        digraph.add_edge('reac%d' % reac_idx, 'spec%d' % spec_idx, weight=filtered_matrix[spec_idx, reac_idx])
                    if filtered_matrix[spec_idx, reac_idx] < 0:
                        digraph.add_edge('spec%d' % spec_idx, 'reac%d' % reac_idx, weight=-filtered_matrix[spec_idx, reac_idx])
        self.tune_graph(digraph, percentage, maximum)


        return digraph      

    def filter_reactions_matrix(self, matrix, maximum):
        # Remove values below threshold
        matrix[np.abs(matrix) < maximum * self.threshold] = 0
        filtered_matrix = np.zeros(np.shape(matrix))
        if self.matrices.analysis in ['RtoS', 'StoRtoS']:
            #Keep reactions involving target species
            for target_spec_idx in self.target_species_idx:
                filtered_matrix[target_spec_idx, :] = matrix[target_spec_idx, :]
        return filtered_matrix


    def filter_species_matrix(self, matrix, maximum):
        '''
        Filtering DIC matrix to remove unnecessary species for graph

        :param DIC: Complete DIC matrix
        :param maxim: Maximum of the DIC over all samples
        :param threshold: Minimum link weight to keep

        :return subDIC: Filtered DIC matrix

        '''
        # Remove values below threshold
        matrix[matrix < maximum * self.threshold] = 0
        filtered_matrix = np.zeros(np.shape(matrix))
        if self.matrices.analysis in ['Glob', 'Pol']:
            digraph = nx.from_numpy_array(matrix, create_using=nx.DiGraph)
            # Get the nodes reachable from the source
            reachable_nodes = [nx.shortest_path(digraph, idx).keys() for idx in self.target_species_idx]
            # Flatten the list
            reachable_nodes = [nodes for sublist in reachable_nodes for nodes in sublist]
            # Keep links only between reachable species
            for node1 in reachable_nodes:
                for node2 in reachable_nodes:
                    filtered_matrix[node1, node2] = matrix[node1, node2]
        elif self.matrices.analysis in ['StoS']:
            # Keep links entering or leaving the target species
            for target_spec_idx in self.target_species_idx:
                filtered_matrix[:,target_spec_idx] = matrix[:,target_spec_idx]
                filtered_matrix[target_spec_idx, :] = matrix[target_spec_idx, :]
        elif self.matrices.analysis in ['StoR', 'StoRtoS']:
            #Keep reactions involving target species
            for target_spec_idx in self.target_species_idx:
                filtered_matrix[target_spec_idx, :] = matrix[:,target_spec_idx]          
        return filtered_matrix

    def tune_graph(self, digraph, percentage, maximum):
        """
        Tune nodes and edges attributes to prepare for plotting
        :param nx_graph_object: Graph to tune (networkx graph object)
        :species: Targeted species

        : return None:
        """
        # Change graph display ratio
        digraph.graph['ratio'] = self.ratio
        # Change nodes properties
        for node, attributes in digraph.nodes(data=True):
            if 'label' in attributes.keys():
                nodename = attributes['label']
            else:
                nodename = node
            if nodename in self.target_species:
                attributes['shape'] = 'doublecircle'
            elif nodename in self.spec_label.values():
                attributes['shape'] = 'circle'
            elif nodename in self.reac_label.values() :
                attributes['shape'] = 'box'
        # Change edges properties
        for node1,node2,attributes in digraph.edges(data=True):
            # Labels and line width
            if percentage:
                attributes['label'] = '%.2f' % attributes.get('weight','') + '%'
                attributes['penwidth'] = 1
            else:
                attributes['label'] = 'w = %.2e' % attributes.get('weight','')
                attributes['penwidth'] = str(abs(attributes.get('weight','')) / maximum * 5)
            # Color
            attributes['color'] = "0.000 1.000 1.000 0.0011578054962364537"


    def dump_graphs(self, digraphs=None, digraphs_integrated=None, directory=None, filename=None):
        if not digraphs:
            digraphs = self.digraphs
        if not digraphs_integrated:
            digraphs_integrated = self.digraphs_integrated
        if not directory:
            directory = os.path.join(analysis_dir, '%s_analysis' % self.matrices.analysis)
        fileio.create_dir(directory, overwrite=False)
        # for case in digraphs:
        #     case_directory = os.path.join(directory, case.nickname)
        #     fileio.create_dir(case_directory, overwrite=False)
        #     for i_sample, grid in enumerate(self.matrices.data_dict[case]['grid'][:]):
        #         self.dump_graph(digraphs[case][i_sample], case_directory, '%s_%.8f' % (self.matrices.analysis, grid), fmt='jpg')
        for case in digraphs_integrated:
            case_directory = os.path.join(directory, case.nickname)
            fileio.create_dir(case_directory, overwrite=False)
            self.dump_graph(digraphs_integrated[case], case_directory, self.matrices.analysis + '_integrated', fmt='jpg')

        return


    def dump_graph(self, digraph, output_folder, file_name, fmt):
        """
        """
        # Write dot file
        nx.drawing.nx_agraph.write_dot(digraph, '%s/%s' % (output_folder, file_name))
        # Read dot file with graphviz
        gv_graph_object = gv.Source.from_file('%s/%s' % (output_folder, file_name))
        if fmt is not None:
            gv_graph_object.render('%s/%s' % (output_folder, file_name), format=fmt)