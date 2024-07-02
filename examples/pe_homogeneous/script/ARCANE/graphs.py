"""Functions for management of validation graphs"""

import csv
import importlib
import math
import os.path
import sys

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

import ARCANE.cases as cases
import ARCANE.display as display
import ARCANE.error as error
import ARCANE.database as database
import ARCANE.kwdict as kwdict
import ARCANE.mechanisms as mechanisms
import ARCANE.tools as tools

logger = display.Logger()
logger.set_log('logCompute')

kwdict = kwdict.Kwdict()


def init_graphs_database(graphdir="graphs", overwrite=False):
    """Initialize graphs folder
    For now, do nothing if it exists. Might include cleanup routines later on.

    Parameters
    ----------
    graphdir : str
        name of folder that will contain all the graphs (Default value = "graphs")
    overwrite : bool
        if True, folder is overwritten (Default value = False)


    Created: 18/04/09 [QC]

    Last modified: 18/04/09 [QC]

    """

    # Initialize case directory
    database.create_dir(graphdir, overwrite)

    # Set up graphdir for all instances
    Graph.graphdir = graphdir


class Graph:

    def __init__(self, style="default"):
        """Main Graph class, containing all class-nonspecific methods and attributes

        Parameters
        ----------
        style: str
            plotting style to be used (Default value = "default")


        Created: 18/04/09 [QC]

        Last modified: 19/08/06 [QC]
        """

        # Default plotting style
        self.style_used = style
        self.plotting_style(style)

        # Init graph directory
        if not hasattr(self, "graphdir"):
            init_graphs_database()
            logger.debug("Initiating default graphs directory : graphs")
        self.image_dir = "images"

        # Predefined colours in html
        self.pretty_colours = {'red': '#A60628',
                               'blue': '#348ABD',
                               'yellow': '#FFC325',
                               'green': '#00b300',
                               'purple': '#9900cc',
                               'salmon': '#ff4d4d',
                               'electric': '#0066ff'}

        # Data extension
        self.extend_data = False

    def show(self):
        """Show the `matplotlib.pyplot` figure """
        plt.show()

    def plotting_style(self, name='rainbow'):
        """Sets the style of plot for following graphs

        Parameters
        ----------
        name : str
            name of the style, either 'paper', 'paper_markers', 'rainbow', 'rainbow_markers', 'default', or path to
            a custom style sheet (Default value = 'rainbow')

        See Also
        --------
        matplotlib.style
            Configuration of `matplotlib.style <https://matplotlib.org/stable/api/style_api.html>`_

        """

        if name not in ['paper', 'paper_markers', 'rainbow', 'rainbow_markers', 'default']:
            plt.style.use(name)
            self.style_used = name

        else:

            path_to_init = tools.__file__

            terminator = path_to_init.rfind('/tools.py')
            path = (path_to_init[:terminator])
            path += '/style_sheets'

            if name == 'paper':

                plt.style.use(path + '/black_paper.mplstyle')
                self.style_used = name

            elif name == 'paper_markers':

                plt.style.use(path + '/black_paper_markers.mplstyle')
                self.style_used = name

            elif name == 'rainbow':

                plt.style.use(path + '/rainbow_paper.mplstyle')
                self.style_used = name

            elif name == 'rainbow_markers':

                plt.style.use(path + '/rainbow_paper_markers.mplstyle')
                self.style_used = name

            else:

                self.style_used = "default"

            logger.debug('Plotting style: ' + self.style_used)

            importlib.reload(plt)

    def arrow(self, x, y, dx, dy, head_width=0.5, head_length=3, show=False, save=None, hold=False):
        """Plot an arrow in the `matplotlib.pyplot` figure according to `arrow module
        <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.arrow.html>`_

        Parameters
        ----------
        x :
            x coordinate
        y :
            y coordinate
        dx :
            x component of the arrow
        dy :
            y component of the arrow
        head_width : float
            head width (Default value = 0.5)
        head_length :
            head length (Default value = 3)
        show : bool
            if True displays the figure (Default value = False)
        save : str
            path to save the file (Default value = None)
        hold : bool
            hold on the graph (Default value = False)

        See Also
        --------
        matplotlib.pyplot.arrow

        """
        plt.arrow(x, y, dx, dy, fc="k", ec="k", head_width=head_width, head_length=head_length)

        if show and not hold:
            plt.show()

        if save:
            plt.savefig(save)

        if not hold:
            plt.clf()

    def quick_graph(self, xdata, ydata, xlabel='', ylabel='', ylim=0, title=None, vertical_line=None, legend=None,
                    marker=None, colour=None, line_style=None, reverse=False, xlog=False, ylog=False, show=True,
                    save=None, hold=False, clip=None, vertical_line_colour=None):
        """Quick plotting with basic `matplotlib` functionalities

        Parameters
        ----------
        xdata :
            abscissa of the graph
        ydata :
            ordinate of the graph
        xlabel :str
            abscissa axis label (Default value = '')
        ylabel : str
            ordinate axis label (Default value = '')
        title : str
            title of the graph (Default value = None)
        vertical_line :
            draws a verticle line at given coordinates (Default value = None)
        legend : str
            legend of the curve (Default value = None)
        symbol : str
            symbol
        colour : str
            colour (Default value = None)
        reverse : bool
            inverts the x axis (Default value = False)
        xlog : bool
            if True x in log scale (Default value = False)
        ylog : bool
            if True y in log scale (Default value = False)
        show : bool
            if True displays the figure (Default value = True)
        save : str
            path to save the file (Default value = None)
        hold : bool
            hold on the graph (Default value = False)
        ylim : float
            maximum value of y axe (Default value = 0). If 0, no ylim is applied
        marker : str
            marker style, typically 'o', 'x', '<', '>'.
            See `matplotlib.markers <https://matplotlib.org/stable/api/markers_api.html>`_ (Default value = None)
        line_style : float
            line style, typically '-', '-.', see possible line style
            `here <https://matplotlib.org/3.5.0/gallery/lines_bars_and_markers/linestyles.html>`_ (Default value = None)
        clip : list
           list of two coordinates :math:`[x_{min}, x_{max}]` to clip the coordinates (Default value = None)
        vertical_line_colour :
            color of vertical line (Default value = None)

        See Also
        --------
        matplotlib.figure

        matplitlib.colors

        matplotlib.markers

        """
        if type(ydata) == list and type(ydata[0]) == list:
            data_list = ydata
        else:
            data_list = [ydata]

        for ydata in data_list:

            ydata = np.array(ydata)

            if not plt.fignum_exists(1):
                plt.figure()

            if clip:
                x_min = min(range(len(xdata)), key=lambda i: abs(xdata[i] - clip[0]))
                x_max = min(range(len(xdata)), key=lambda i: abs(xdata[i] - clip[1]))
            else:
                x_min = 0
                x_max = -1

            xdata = xdata[x_min:x_max]
            ydata = ydata[x_min:x_max]

            plt.plot(xdata, ydata, label=legend, color=colour, marker=marker, linestyle=line_style)

        if reverse and xdata[0] != xdata[-1]:
            plt.xlim(xdata[0], xdata[-1])

        if vertical_line:
            if not isinstance(vertical_line, list):
                vertical_line = [vertical_line]
            for value in vertical_line:
                plt.axvline(x=value, color=vertical_line_colour)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if title:
            plt.title(title)

        if legend:
            plt.legend()

        if xlog:
            plt.xscale('log')

        if ylog:
            plt.yscale('log')

        if ylim != 0:
            plt.ylim(0, ylim)

        if save:
            plt.savefig(save)
            logger.info('Figure saved as ' + save)

        if show:
            plt.show()

        if not hold:
            plt.clf()

    def case_plot(self, cases_list, mechanisms_list, x, y, normalized=False, offset=False,
                  xlabel=None, ylabel=None,
                  show=True, save=False, hold=False, legend=None, title=None, num_fig=1,
                  clip=None, mark_every=1, xlog=False, ylog=False):
            """
            Quick plotting of profiles

            Parameters
            ----------
            cases_list : list
                list of class :func:`~ARCANE.cases.Case` objects
            mechanisms_list : list
                list of class :func:`~ARCANE.mechanisms.Mechanism` objects
            x : str
                name of the variable wanted as an abscissa
            y : str
                name of the variable wanted as an ordinate or list of names
            normalized : bool or str
                if True, the profiles will be normalized by their maximum (Default value = False)
                if a string corresponding to a method (max, min, mean, etc.) is provided it will serve for the normalization
            offset : bool or str
                if True, the profiles will be subtracted by their mean (Default value = False)
                if a string corresponding to a method (max, min, mean, etc.) is provided it will serve as offset
            show : bool
                if False, doesn't display the figure (Default value = True)
            save : bool
                if True, saves the figure in pdf format (Default value = False)
            hold : bool
                keep the previous graphs in memory (Default value = False)
            legend : str
                gives a non default label to the plot
            title : str
                custom title for the graph (Default value = None)
            num_fig : int
                number of the figure (Default value = 1)
            clip : list
                list of two coordinates [x_min, x_max] to clip the coordinates (Default value = None)
            xlabel : str
                abscissa axis label (Default value = None)
            ylabel : str
                ordinate axis label (Default value = None)
            mark_every :
                `markevery` argument (Default value = 1)
            xlog : bool
                if True x in log scale (Default value = False)
            ylog : bool
                if True x in log scale (Default value = False)

            Created: 18/04/09 [QC]

            Last modified: 18/04/09 [QC]

            """
            if not cases_list:
                logger.error("No case was specified")
                return

            if not mechanisms_list:
                logger.error("No mechanism was specified")
                return

            if not isinstance(cases_list, list):
                cases_list = [cases_list]

            if not isinstance(mechanisms_list, list):
                mechanisms_list = [mechanisms_list]

            if title or title is False:
                auto_title = False
            else:
                auto_title = True

            plt.figure(num_fig)

            if x in kwdict.names['grid']:
                x = kwdict.reactor_grid_name[cases_list[0].reactor_type]

            for case in cases_list:

                for index_mech, mechanism in enumerate(mechanisms_list):

                    data_dict = case.data_dict(mechanism)

                    # Updating kwdict with the species_names
                    kwdict.update('names', mechanism.kwdict_names_extension)
                    kwdict.update('units', mechanism.kwdict_units_extension)

                    if type(y) == list:
                        test = all(y) not in kwdict.names and all(kwdict.get_base_name(y)) not in data_dict
                    else:
                        test = y not in kwdict.names and kwdict.get_base_name(y) not in data_dict

                    if test:
                        # Updating kwdict with the species_names
                        kwdict.update('names', mechanism.kwdict_names_extension)
                        kwdict.update('units', mechanism.kwdict_units_extension)

                        data_dict = case.data_dict(mechanism, extend_data=True)

                    # Checking if the keyword is correct
                    if x not in data_dict and kwdict.get_base_name(x) not in data_dict:

                        data_dict = case.data_dict(mechanism, reload=True)

                    if x not in data_dict and kwdict.get_base_name(x) not in data_dict:
                        logger.error('The keyword for the abscissa is not a valid one: ' + x)
                        logger.error('Here is a list of the ones you can use: ' + ', '.join(data_dict.keys()))
                        return

                    abscissa = data_dict[kwdict.get_base_name(x)]

                    if clip:
                        x_min = min(range(len(abscissa)), key=lambda i: abs(abscissa[i] - clip[0]))
                        x_max = min(range(len(abscissa)), key=lambda i: abs(abscissa[i] - clip[1]))
                    else:
                        x_min = 0
                        x_max = -1

                    abscissa = abscissa[x_min:x_max]

                    one_warning = False
                    if type(y) == list:

                        for ordi in y:
                            # Checking if the keyword is correct
                            if ordi not in data_dict and kwdict.get_base_name(ordi) not in data_dict:
                                logger.error('List of keywords for the ordinate values are only available for species')
                                logger.error(ordi + ' was given')
                                logger.error(
                                        'Here is a list of the species you can use: ' + ', '.join(
                                                mechanism.ctmech.species_names))
                                return

                            ordinate = data_dict[kwdict.get_base_name(ordi)]
                            ordinate = ordinate[x_min:x_max]

                            if normalized:
                                if isinstance(normalized, str):
                                    normalizer = tools.extract_scalar(ordinate, normalized)

                                elif isinstance(normalized, float):
                                    normalizer = normalized

                                else:
                                    normalizer = tools.extract_scalar(ordinate, 'max')

                                if normalizer > 1e-60:
                                    ordinate = ordinate / normalizer
                                else:
                                    continue

                            if offset:
                                if isinstance(offset, str):
                                    offseter = tools.extract_scalar(ordinate, offset)

                                elif isinstance(normalized, float):
                                    offseter = offset

                                else:
                                    offseter = tools.extract_scalar(ordinate, 'mean')

                                if offseter > 1e-60:
                                    ordinate = ordinate - offseter
                                else:
                                    continue

                            if all(y) in mechanism.species_names:
                                ordinate_name = kwdict.get_full_name(ordi).split(' ')
                                label_name = ordinate_name[0]
                                ordinate_name = ' '.join(ordinate_name[1:])
                            else:
                                label_name = ordi
                                ordinate_name = kwdict.get_full_name(ordi)

                            if len(y) < 20:
                                if not legend and legend != False:
                                    label = label_name + ' '
                                    if len(cases_list) > 1:
                                        label += case.nickname + ' '
                                    if len(mechanisms_list) > 1:
                                        label += mechanism.nickname
                                else:
                                    label = legend
                            elif not one_warning:
                                label = None
                                logger.info('WARNING No legend will be displayed as the number of curves is to high')
                                one_warning = True

                            # Plotting
                            if not label:
                                label = None

                            plt.plot(abscissa, ordinate, label=label)

                            if not ylabel:
                                if not normalized:
                                    plt.ylabel(ordinate_name + ' [' + kwdict.units[ordinate_name] + ']')
                                else:
                                    plt.ylabel('Normalized ' + ordinate_name  + ' [-]')

                                if auto_title:
                                    title = kwdict.reactor_labels[case.reactor_type][0] + ' Mass fraction vs ' \
                                            + ' ' + kwdict.get_full_name(x) + ' at '
                            else:
                                plt.ylabel('ylabel')

                    else:
                        # Checking if the keyword is correct
                        if y not in data_dict and kwdict.get_base_name(y) not in data_dict:
                            # Updating kwdict with the species_names
                            kwdict.update('names', mechanism.kwdict_names_extension)
                            kwdict.update('units', mechanism.kwdict_units_extension)

                            data_dict = case.data_dict(mechanism, extend_data=True)

                        if y not in data_dict and kwdict.get_base_name(y) not in data_dict:
                            logger.error('The keyword for the ordinate is not a valid one: ' + y)
                            logger.error('Here is a list of the ones you can use: ' + ', '.join(data_dict.keys()))
                            return

                        ordinate = data_dict[kwdict.get_base_name(y)]
                        ordinate = ordinate[x_min:x_max]

                        if normalized:
                            if isinstance(normalized, str):
                                normalizer = tools.extract_scalar(ordinate, normalized)

                            elif isinstance(normalized, float):
                                normalizer = normalized

                            else:
                                normalizer = tools.extract_scalar(ordinate, 'max')

                            if normalizer > 1e-60:
                                ordinate = ordinate / normalizer
                            else:
                                continue

                        if offset:
                            if isinstance(offset, str):
                                offseter = tools.extract_scalar(ordinate, offset)

                            elif isinstance(normalized, float):
                                offseter = offset

                            else:
                                offseter = tools.extract_scalar(ordinate, 'mean')

                            if offseter > 1e-60:
                                ordinate = ordinate - offseter
                            else:
                                continue

                        if not legend and legend is not False:
                            label = ''
                            if len(cases_list) > 1:
                                label += case.nickname + ' '
                            if len(mechanisms_list) > 1:
                                label += mechanism.nickname
                        else:
                            label = legend

                        # Plotting
                        plt.plot(abscissa, ordinate, label=label,
                                 color=mechanism.colour, marker=mechanism.marker, linestyle=mechanism.line_style,
                                 markevery=mark_every)

                        # plt.axvline(x=40e-6, color='grey', alpha=0.5)

                        if not ylabel:
                            y = y.replace('_ptcl', '')
                            if not normalized:
                                plt.ylabel(kwdict.get_full_name(y) + ' [' + kwdict.units[kwdict.get_full_name(y)] + ']')
                            else:
                                plt.ylabel('Normalized ' + kwdict.get_full_name(y) + ' [-]')
                        else:
                            plt.ylabel(ylabel)

                        if auto_title:
                            title = kwdict.reactor_labels[case.reactor_type][0] + ' at '

                if not xlabel:
                    x = x.replace('_ptcl', '')
                    plt.xlabel(kwdict.get_full_name(x) + ' [' + kwdict.units[kwdict.get_full_name(x)] + ']')
                else:
                    plt.xlabel(xlabel)

            if auto_title and len(cases_list) == 1:
                if hasattr(case, 'temperature'):
                    title += 'T = ' + str(case.temperature) + ' K '
                if hasattr(case, 'pressure'):
                    title += 'P = ' + str(case.pressure / 1e5) + ' bar '
                if hasattr(case, 'phi'):
                    title += 'phi = ' + str(case.phi)

            if title:
                plt.title(title)

            if xlog:
                plt.xscale('log')

            if ylog:
                plt.yscale('log')

            if label:
                plt.legend()

            if ylog:
                plt.yscale('log')

            if type(y) == list:
                y = 'Y'

            if save and not hold:

                if type(save) == str:
                    plt.savefig(save)
                    logger.info('Graph saved as ' + save)

                else:
                    plt.savefig(self.graphdir + '/' + 'X' + x + 'Y' + y + mechanism.nickname + '.pdf')
                    logger.info('Graph saved as ' + self.graphdir + '/' + 'X' + x + 'Y' + y + mechanism.nickname + '.pdf')

            if show and not hold:
                plt.show()

            if not hold:
                plt.clf()

    def parametric_plot(self, cases_list, mechanisms_list, x, y, reactor_type=None, P=-1, T=-1, phi=-1, fuel={},
                        xlog=None, ylog=None, show=True, save=False, hold=False, legend=None, title=None,
                        error_plot=None, num_fig=1, normalized=False, xlabel=None, ylabel=None):
        """Plots a parametric plot for 2 fixed parameters between P, T and phi

        Parameters
        ----------
        cases_list : list
            list of class :func:`~ARCANE.cases.Case` object
        mechanisms_list : list
            list of class :func:`~ARCANE.mechanisms.Mechanism` object
        x : str
            name of the variable wanted as an abscissa
        y : str
            name of the variable wanted as an ordinate
        reactor_type :
            type of reactor (Default value = None)
        P : float
            fixed pressure parameter (Default value = -1)
        T : float
            fixed temperature parameter (Default value = -1)
        phi : float
            fixed equivalence parameter (Default value = -1)
        fuel : dict
            fuel composition parameter (Default value = {})
        show : bool
            if False, doesn't display the figure (Default value = True)
        xlog : bool
            if True, logarithmic scale is applied on abscissa (Default value = None)
        ylog : bool
            if True, logarithmic scale is applied on ordinate (Default value = None)
        save : bool
            if True, saves the figure in pdf format (Default value = False)
        hold : bool
            keep the previous graphs in memory (Default value = False)
        legend : str
            gives a non default label to the plot
        title : str
            custom title for the graph (Default value = None)
        error_plot : bool
            if True, the mechanisms of the list are compared to the first one (Default value = None)
        num_fig : int
            number of the figure (Default value = 1)
        normalized : str | float | bool
            the profiles can be normalized using the `ARCANE.kwdict.methods` or a given value. If True, the maximum is
            used(Default value = False)
        xlabel : str
            abscissa axis label (Default value = None)
        ylabel : str
            ordinate axis label (Default value = None)

        Returns
        -------
        x_data : list
            data of x-axis
        y_axis : list
            data of y-axis

        """

        # Preventing the non backward compatibility with former mandatory arguments of the function
        if not isinstance(mechanisms_list, list):
            if error_plot:
                logger.error('Error! A list of mechanisms is required in order to compute the error.')
                return

            mechanisms_list = [mechanisms_list]

        plt.figure(num_fig)

        if reactor_type or P > -1 or T > -1 or phi > -1 or fuel:
            selected_cases = tools.get_case_by_state(cases_list, reactor_type=reactor_type, P=P, T=T, phi=phi, fuel=fuel)
        else:
            selected_cases = cases_list

        if len(selected_cases) == 0:
            logger.warning('WARNING: No matching cases for parametric plot')
            return

        base_value = False
        if x in kwdict.names['T'] + kwdict.names['1000/T'] + kwdict.names['P'] + kwdict.names['phi'] \
            or x.startswith('Fuel'):
            base_value = True
            x_quantity = x
            x_method = ''
        else:
            # Extracting quantity and method from string for x
            x_split = x.split(' ')
            x_last_split = x_split[-1]
            x_quantity = x
            x_method = ''
            # Checking if last word is a method
            for name in kwdict.methods:
                if x_last_split in kwdict.methods[name]:
                    x_method = x_last_split
                    x_quantity = ' '.join(x_split[:-1])

            # Checking if last word is a compute method
            for name in kwdict.compute_methods:
                if x_last_split in kwdict.compute_methods[name]:
                    x_method = x_last_split
                    x_quantity = ' '.join(x_split[:-1])

        # Extracting quantity and method from string for y
        y_split = y.split(' ')
        y_last_split = y_split[-1]
        y_quantity = y
        y_method = ''
        # Checking if last word is a method
        for name in kwdict.methods:
            if y_last_split in kwdict.methods[name]:
                y_method = y_last_split
                y_quantity = ' '.join(y_split[:-1])

        # Checking if last word is a compute method
        for name in kwdict.compute_methods:
            if y_last_split in kwdict.compute_methods[name]:
                y_method = y_last_split
                y_quantity = ' '.join(y_split[:-1])

        if error_plot:
            iter_mech = mechanisms_list[1:]
        else:
            iter_mech = mechanisms_list

        x_data = []
        y_data = []

        verification = True

        # Plotting loop
        for index_mech, mechanism in enumerate(iter_mech):
            abscissa = []
            ordinate = []

            kwdict.update('names', mechanism.kwdict_names_extension)
            kwdict.update('units', mechanism.kwdict_units_extension)

            last_valid_case = None

            for case in selected_cases:

                data_dict = case.data_dict(mechanism)
                test = y not in kwdict.names and kwdict.get_base_name(y) not in data_dict \
                       and kwdict.get_base_name(y) not in ['sl', 'tig', 'thickness']
                if test:
                    # Updating kwdict with the species_names
                    kwdict.update('names', mechanism.kwdict_names_extension)
                    kwdict.update('units', mechanism.kwdict_units_extension)
                    data_dict = case.data_dict(mechanism, extend_data=True)

                case.axis = []

                if not base_value:
                    abscissa.append(case.extract_quantity(x, mechanism=mechanism))
                else:
                    if x in kwdict.names['T']:
                        abscissa.append(case.temperature)

                    elif x in kwdict.names['1000/T']:
                        abscissa.append(1000 / case.temperature)

                    elif x in kwdict.names['P']:
                        abscissa.append(case.pressure)

                    elif x in kwdict.names['phi']:
                        abscissa.append(case.phi)

                    elif x.startswith('Fuel'):
                        logger.error("Not implemented yet ...")

                    else:
                        logger.error("You should not end up here ... That's awkward !")

                if database.check_solution_existence(case, mechanism):

                    if verification:
                        data_names = list(data_dict.keys())
                        data_names += kwdict.names['sl'] + kwdict.names['tig'] + kwdict.names['thickness']

                        # Checking if the keyword is correct
                        if x_quantity not in data_names and kwdict.get_base_name(
                                x_quantity) not in data_names and not base_value:
                            logger.error('The keyword for the abscissa quantity is not a valid one: ' + x_quantity)
                            logger.error('Here is a list of the ones you can use: ' + ', '.join(data_names))
                            return

                        # Checking if the keyword is correct
                        if y_quantity not in data_names and kwdict.get_base_name(y_quantity) not in data_names:
                            logger.error('The keyword for the ordinate quantity is not a valid one: ' + y_quantity)
                            logger.error('Here is a list of the ones you can use: ' + ', '.join(data_names))
                            return

                        # Checking if the keyword is correct
                        possible_methods = ['']
                        for name in kwdict.methods:
                            possible_methods += kwdict.methods[name]

                        for name in kwdict.compute_methods:
                            possible_methods += kwdict.compute_methods[name]

                        if x_method not in possible_methods and not base_value:
                            logger.error('The keyword for the abscissa method is not a valid one: ' + x_method)
                            logger.error('Here is a list of the ones you can use: ' + ', '.join(possible_methods))
                            return

                        if y_method not in possible_methods:
                            logger.error('The keyword for the ordinate method is not a valid one: ' + y_method)
                            logger.error('Here is a list of the ones you can use: ' + ', '.join(possible_methods))
                            return

                        verification = False

                    if error_plot:
                        error_value, ref, curr = error.compute_error(y, case, mechanisms_list[0], mechanism=mechanism)
                        if error_value * 100 <= 100:
                            ordinate.append(error_value * 100)
                        else:
                            ordinate.append(np.nan)
                    else:
                        ordinate.append(case.extract_quantity(y, mechanism=mechanism))
                else:
                    logger.debug(f"{case.solution} does not exists so its value will be NaN")
                    ordinate.append(np.nan)

            if not legend and legend is not False:
                label = mechanism.nickname
            else:
                label = legend

            if normalized:
                if isinstance(normalized, str):
                    normalizer = tools.extract_scalar(ordinate, normalized)

                elif isinstance(normalized, float):
                    normalizer = normalized

                else:
                    normalizer = tools.extract_scalar(ordinate, 'max')

                if normalizer > 1e-60:
                    ordinate = ordinate / normalizer
                else:
                    continue

            plt.plot(abscissa, ordinate, label=label,
                     color=mechanism.colour, marker=mechanism.marker, linestyle=mechanism.line_style)

            x_data.append(abscissa)
            y_data.append(ordinate)

        if y in kwdict.names['tig'] and ylog == None:
            ylog = True
        if y in kwdict.names['SootVolumeFraction'] and ylog == None:
            ylog = True
        if y in kwdict.names['SootNumberDensity'] and ylog == None:
            ylog = True

        if xlog:
            plt.xscale('log')

        if ylog:
            plt.yscale('log')

        if not xlabel:
            if x_quantity in mechanism.ctmech.species_names:
                if error_plot:
                    xlabel = 'Error on ' + kwdict.get_full_name(x_method) + ' ' + x_quantity + ' mass fraction ' + ' [%]'
                else:
                    xlabel = kwdict.get_full_name(x_method) + ' ' + x_quantity + ' mass fraction ' + ' [-]'
            else:
                if error_plot:
                    xlabel = 'Error on ' + kwdict.get_full_name(x_method) + ' ' + kwdict.get_full_name(x_quantity) + ' [%]'
                else:
                    xlabel = kwdict.get_full_name(x_method) + ' ' + kwdict.get_full_name(x_quantity) \
                             + ' [' + kwdict.units[kwdict.get_full_name(x_quantity)] + ']'

        if not ylabel:
            if y_quantity in mechanism.ctmech.species_names:
                if error_plot:
                    ylabel = 'Error on ' + kwdict.get_full_name(y_method) + ' ' + y_quantity + ' mass fraction ' + ' [%]'
                else:
                    ylabel = kwdict.get_full_name(y_method) + ' ' + y_quantity + ' mass fraction ' + ' [-]'
            else:
                if error_plot:
                    ylabel = 'Error on ' + kwdict.get_full_name(y_method) + ' ' + kwdict.get_full_name(y_quantity) + ' [%]'
                else:
                    ylabel = kwdict.get_full_name(y_method) + ' ' + kwdict.get_full_name(y_quantity) \
                             + ' [' + kwdict.units[kwdict.get_full_name(y_quantity)] + ']'

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if title:
            plt.title(title)

        if label:
            plt.legend()

        if save:  # and not hold:

            if type(save) == str:

                plt.savefig(save, dpi=2000)
                logger.info('Figure saved to ' + save)

            else:

                path_to_save = self.graphdir + '/Parametric' + 'X' + x.replace(' ', '_').replace('/', '-') \
                               + 'Y' + y.replace(' ', '_') + '.pdf'

                if reactor_type:
                    path_to_save += reactor_type
                if P:
                    path_to_save += 'P' + str(P / 1e5)
                if T:
                    path_to_save += 'T' + str(T)
                if phi:
                    path_to_save += 'phi' + str(phi)
                if fuel:
                    path_to_save += 'fuel' + str(fuel)

                path_to_save += '.pdf'

                plt.savefig(path_to_save)
                logger.info('Figure saved to ' + path_to_save)

        if show and not hold:
            plt.show()

        if not hold:
            plt.clf()

        return x_data, y_data

    def plot_timescales(self, cases_list, mechanisms_list, time_step_exponent=0, method='explicit',
                        abs_tol=1e-9, rel_tol=1e-5, species=None,
                        show=True, save=False):
        r"""Plots the timescales of species in a bar plot

        Parameters
        ----------
        cases_list : list
            list of class :func:`~ARCANE.cases.Case` object objects
        mechanisms_list : list
            list of class :func:`~ARCANE.mechanisms.Mechanism` object object
        time_step_exponent : int
            simulation timescale exponent to display a line on the graph (Default value = 0)
        method : str
            - 'explicit': :math:`c / \dot{c}`, useful for explicit codes as it gives an idea of the time step necessary for a proper tempora resolution,

            - 'implicit' :math:`\mathrm{d}c / \mathrm{d}\dot{c}`, represent the real time scale of the species and is the one used in kinetics analysis, in the selection of QSS species for example,

            - 'both', will apply both 'explicit' and 'implicit' methods
            (Default value = 'explicit')
        abs_tol : float
            absolute tolerance, if the concentration of the species is inferior to abs_tol
            its timescale will not be computed (Default value = 1e-9)
        rel_tol : float
            relative tolerance (Default value = 1e-5)
        species : list
            list of the specific species to plot (Default value = None)
        show : bool
            if False, doesn't display the figure (Default value = True)
        save : bool
            if True, saves the figure in pdf format (Default value = False)

        """

        if type(cases_list) != list:
            cases_list = [cases_list]

        if type(mechanisms_list) != list:
            mechanisms_list = [mechanisms_list]

        if species and type(species) != list:
            species = [species]

        n_mech = len(mechanisms_list)

        min_timescale_log = []

        # Telling the user which method is used to compute the timescales
        if method == 'both':
            plural = 's'
        else:
            plural = ''

        logger.info(f"The time scales are computed with {method} method{plural}.")

        # Pre-processing the chosen species
        if not species:
            all_species = []
            for mech in mechanisms_list:
                all_species += mech.species_names
            all_species = list(set(all_species))
        else:
            all_species = species

        plt.figure(figsize=(20, 10))

        for index, mechanism in enumerate(mechanisms_list):

            if method == 'implicit' or method == 'both':
                timescales_dict = tools.get_species_timescale(cases_list, mechanism, abs_tol=abs_tol, rel_tol=rel_tol)

                # species_names = timescales_dict.keys()
                # min_timescales = timescales_dict.values()
                min_timescales = [timescales_dict[spec] if spec in timescales_dict else 0 for spec in all_species]

                min_timescale_log_im = [- np.log10(x) if x > 0 else 0 for x in min_timescales]

                min_timescale_log += min_timescale_log_im

                x = range(len(all_species))
                ns = len(all_species)

                if method == 'both':
                    colour = 'k'
                    width = 0.8 / n_mech
                else:
                    colour = mechanism.colour
                    width = 0.8 / n_mech

                x = np.array(x)
                x = x - 0.4 + ((1/n_mech)*0.4) * (2 * index + 1)

                label = None
                if n_mech > 1:
                    label = f"{mechanism.nickname}"
                    if method == 'both':
                        label += ' implicit'
                else:
                    if method == 'both':
                        label = 'implicit'

                plt.bar(x, min_timescale_log_im, width, color=colour, label=label)

            if method == 'explicit' or method == 'both':
                timescales_dict = tools.get_species_timescale(cases_list, mechanism, explicit=True,
                                                              abs_tol=abs_tol, rel_tol=rel_tol)

                # species_names = timescales_dict.keys()
                # min_timescales = timescales_dict.values()
                min_timescales = [timescales_dict[spec] if spec in timescales_dict else 0 for spec in all_species]

                min_timescale_log_ex = [- np.log10(x) if x > 0 else 0 for x in min_timescales]

                min_timescale_log += min_timescale_log_ex

                x = range(len(all_species))
                ns = len(all_species)
                width = 0.8 / (n_mech*1.5)
                x = np.array(x)
                x = x - 0.4 + ((1/n_mech)*0.4) * (2 * index + 1)

                label = None
                if n_mech > 1:
                    label = f"{mechanism.nickname}"
                    if method == 'both':
                        label += ' explicit'
                else:
                    if method == 'both':
                        label = 'explicit'

                plt.bar(x, min_timescale_log_ex, width, color=mechanism.colour, label=label)

            if method not in ['both', 'explicit', 'implicit']:
                logger.error("ERROR ! Wrong keyword for the method, try 'implicit', 'explicit' or 'both'")
                return

        time_step_exponent = int(abs(time_step_exponent))

        max_exponent = int(max(min_timescale_log))

        exponents = list(range(time_step_exponent, max_exponent + 2))

        if time_step_exponent != 0:
            for expo in exponents:
                plt.plot(x, expo * np.ones(ns), 'r--')

        x_abs = np.array(range(len(all_species)))
        plt.xticks(x_abs, all_species, rotation=90)
        plt.xlabel('Species of interest')
        plt.ylabel(r'$-\log(\tau_c)$')
        plt.legend()
        plt.tight_layout()

        name_for_saving = mechanism.nickname.replace('/', '_')

        if save:
            if type(save) == str:
                plt.savefig(save)
                logger.info('Graph saved as ' + save)
            else:
                plt.savefig(self.graphdir + '/' + 'Timescales_' + name_for_saving + '.pdf')
                logger.info('Graph saved as ' + self.graphdir + '/' + 'Timescales_' + name_for_saving + '.pdf')

        if show:
            plt.show()

    def ternary_contour_from_cases_list(self, cases_list, mechanism, reactor_type, value_name, list_blend=[],
            P=-1, T=-1, phi=-1, save_image=True, show=False, error_plot=False, reference_mechanism=None):
        """Function to build a ternary contour plot from a cases_list for a given variable, requires a 3 species surrogate.

        Parameters
        ----------
        cases_list : list
            list of class :func:`~ARCANE.cases.Case` objects in which the solution are defined
        mechanism :  class :func:`~ARCANE.mechanisms.Mechanism` object
            mechanism object from class :func:`~ARCANE.mechanisms.Mechanism` object class
        reactor_type : str
            type of reactor, `ARCANE.kwdict.reactor_labels`, typically 'C0DV', 'C0DP', 'C1DP', 'C1DCFD', 'C1DF'
        value_name : str
            name of the variable wanted as the contour field (from error dict)
        list_blend : list
            list of the 3 species in the fuel blend. if not given, the function finds them itself (Default value = [])
        P : float
            fixed pressure parameter (Default value = -1)
        T : float
            fixed temperature parameter (Default value = -1)
        phi : float
            fixed equivalence parameter (Default value = -1)
        save_image : bool
            if True, create a sub-repertory images in the current folder if it does not already exist, and save figure in it (Default value = True)
        show : bool
            if False, doesn't display the figure (Default value = False)
        error_plot : bool
            if True, compute the error compared to a reference mechanism (Default value = False)
        reference_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            reference mechanism required to plot an error_plot (Default value = None)


        See Also
        --------
        plotly.figure_factory


        Created: 23/03/2021 [TL]

        """
        import os.path

        import plotly.figure_factory as ff

        # Check if the cases_list given contains only fuel with 3 different species
        ##########################################################################
        # Recup the fuel species for all the cases
        if not list_blend:
            list_blend = []
            fuel_species = []
            for case in cases_list:
                if str(list(case.fuel.keys())) not in fuel_species:
                    fuel_species.append(str(list(case.fuel.keys())))

            fuel_species = [eval(fuel_string) for fuel_string in fuel_species]

            # From the fuel species of all the cases, extract only species that are different

            for case in fuel_species:
                for species in case:
                    if species not in list_blend:
                        list_blend.append(species)

        # Check if there is only three fuel species in the cases_list given
        if len(list_blend) != 3:
            logger.error("The cases_list given does not contain a 3 species fuel blend")
            return

        # Check if thermo data are given
        if P == -1 or T == -1 or phi == -1:
            logger.error("Specified (T,P,phi) must be given to compute a ternary contour plot")
            return

        # If error_plot is true, then a reference mechanism must be provided
        if error_plot:
            if not reference_mechanism:
                logger.error("For an error plot, a reference mechanism must be provided")
                return

        # Extract the data from the cases_list
        ####################################
        # Get the cases for the specified reactor and thermo conditions in the cases_list
        selected_cases = tools.get_case_by_state(cases_list, reactor_type, P=P, T=T, phi=phi)

        # Initialize list to get the values of interest
        X_fuel_a = []
        X_fuel_b = []
        X_fuel_c = []
        value = []

        # Go through the cases selected
        for case in selected_cases:

            # Initialize check to add a zero molar fraction if species is not in the blend
            check_a = False
            check_b = False
            check_c = False

            # Recup molar fractions from fuel species
            for component in case.fuel.keys():
                if component == list_blend[0]:
                    X_fuel_a.append(case.fuel[component])
                    check_a = True
                if component == list_blend[1]:
                    X_fuel_b.append(case.fuel[component])
                    check_b = True
                if component == list_blend[2]:
                    X_fuel_c.append(case.fuel[component])
                    check_c = True

            # If species is absent from blend, add a zero molar fraction
            if not check_a:
                X_fuel_a.append(0)
            if not check_b:
                X_fuel_b.append(0)
            if not check_c:
                X_fuel_c.append(0)

            # Recup values of interest
            if error_plot:
                error_value, rf, curr = error.compute_error(value_name, case, reference_mechanism, mechanism=mechanism)
                value.append(error_value)
            else:
                value.append(case.extract_quantity(value_name, mechanism=mechanism))

        # Save image
        fig = self.plot_ternary_graph(X_fuel_a, X_fuel_b, X_fuel_c, value, pole_labels=list_blend)

        if save_image:
            if type(save_image) == str:
                fig.write_image(save_image)
                logger.info('Image saved as ' + save_image)

            else:

                if not os.path.exists(self.image_dir):
                    os.mkdir(self.image_dir)

                # Extracting quantity and method from string
                y_split = value_name.split(' ')
                last_split = y_split[-1]

                quantity = value_name
                method = ''

                # Checking if last word is a method
                for name in kwdict.methods:
                    if last_split in kwdict.methods[name]:
                        method = last_split
                        quantity = ' '.join(y_split[:-1])

                value_string = kwdict.get_full_name(quantity) + kwdict.get_full_name(method)
                value_string = value_string.replace(' ', '')

                image_name = f"{self.image_dir}/{value_string}_T{T}_P{P / 1e5}_phi{phi}.pdf"

                fig.write_image(image_name)
                logger.info('Image saved as ' + image_name)

        if show:
            fig.show()

    def ternary_contour_automatic(self, mechanism, reactor_type, value_name, fuel_spec_a, fuel_spec_b, fuel_spec_c,
            T=-1, P=-1, phi=-1, nb_pts=21, show=False, save_image=True, save_data=False, error_plot=False,
            reference_mechanism=None):
        """Function to build a ternary graph for a given variable, requires a 3 species surrogate.
        The function chooses itself the point to plot from a number of points given by the user.

        Parameters
        ----------
        mechanism :
            mechanism object from class :func:`~ARCANE.mechanisms.Mechanism` object
        reactor_type : str
            type of reactor, `ARCANE.kwdict.reactor_labels`, typically 'C0DV', 'C0DP', 'C1DP', 'C1DCFD', 'C1DF'
        value_name : str
            name of the variable wanted as the contour field (from error dict)
        fuel_spec_a :
            name of the first species on the 3 surrogate fuel blend as a string
        fuel_spec_b :
            name of the second species on the 3 surrogate fuel blend as a string
        fuel_spec_c :
            name of the third species on the 3 surrogate fuel blend as a string
        P : float
            fixed pressure parameter (Default value = -1)
        T : float
            fixed temperature parameter (Default value = -1)
        phi : float
            fixed equivalence parameter (Default value = -1)
        nb_pts : int
            number of points to compute for the contour field (choose at least 10 points, 21 by default)
        show : bool
            if False, doesn't display the figure (Default value = False)
        save_image : bool
            if True, create a sub-repertory images in the current folder if it does not already exist, and save figure in it (Default value = True)
        save_data : bool
            if True, write a csv file containing the molar concentrations and the contour field variable (Default value = False)
        error_plot : bool
            if True, compute the error compared to a reference mechanism (Default value = False)
        reference_mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            reference mechanism required to plot an error_plot (Default value = None)


        Created: 22/03/2021 [TL]

        """
        import plotly.figure_factory as ff

        # Check of T,p,phi or given
        if T == -1 or P == -1 or phi == -1:
            logger.error("You should provide specified T, P and phi to compute the cases for ternary plot")
            return

        # Check if fuel species given are in the mechanism
        nb_fuel_spec_in_mech = 0
        for species in mechanism.species_names:
            if species == fuel_spec_a or species == fuel_spec_b or species == fuel_spec_c:
                nb_fuel_spec_in_mech += 1
        if nb_fuel_spec_in_mech < 3:
            logger.error("At least one of the fuel species given is not in the mechanism")
            return

        # If error_plot is True, a reference mechanism must be given
        if error_plot:
            if not isinstance(reference_mechanism, mechanisms.Mechanism):
                logger.error("For an error_plot, the reference mechanism object must be given")
                return

        # From the number of points given, calculate the molar concentration list
        division = int((-1 + math.sqrt(1 + 4 * 2 * nb_pts)) / 2)
        list_fuel_a = np.linspace(0, 1, division)
        list_fuel_b = np.linspace(0, 1, division)

        # Air composition
        air = "X/O2/0.21/N2/0.79"

        # Lists to recup values of interest
        X_fuel_a = []
        X_fuel_b = []
        X_fuel_c = []
        value = []

        # Loop on the different blend molar composition
        for x in list_fuel_a:
            for y in list_fuel_b:
                if (x + y) <= 1:
                    z = abs(1 - x - y)

                    # Define fuel composition
                    fuel = 'X/' + fuel_spec_a + '/' + str(x) + '/' + fuel_spec_b + '/' + str(
                            y) + '/' + fuel_spec_c + '/' + str(z)

                    # Define cases to run
                    cases_list = []
                    cases_list.extend(cases.create_case(reactor=reactor_type,
                                                        mechanism=mechanism,
                                                        fuel=fuel,
                                                        oxidizer=air,
                                                        pressure=str(P),
                                                        temperature=str(T),
                                                        phi=str(phi)))

                    # Run case
                    try:
                        cases.run_cases(cases_list, mechanism, overwrite=False)
                        if error_plot:
                            cases.run_cases(cases_list, reference_mechanism, overwrite=False)
                    except:
                        logger.warning("A case has not converged and will be ignored for the ternary plot")
                        break

                    # Calculate error for error_plot
                    if error_plot:
                        error_value, ref, curr = error.compute_error(value_name, cases_list[0], reference_mechanism,
                                                                     mechanism=mechanism)

                    # Recup defined molar concentration
                    X_fuel_a.append(x)
                    X_fuel_b.append(y)
                    X_fuel_c.append(z)

                    # Recup value of interest for the case
                    if error_plot:
                        value.append(error_value)
                    else:
                        value.append(cases_list[0].extract_quantity(value_name, mechanism=mechanism))

        # Save images
        fig = self.plot_ternary_graph(X_fuel_a, X_fuel_b, X_fuel_c, value,
                                 pole_labels=[fuel_spec_a, fuel_spec_b, fuel_spec_c])

        if save_image:
            if not os.path.exists("images"):
                os.mkdir("images")
            fig.write_image("images/" + value_name + "_T_" + str(T) + "_P_" + str(P) + "_phi_" + str(phi) + ".pdf")

        if show:
            fig.show()

        # Save results
        if save_data:
            fieldnames = [fuel_spec_a, fuel_spec_b, fuel_spec_c, value_name]
            with open('images/data_' + str(mechanism.name) + '.csv', 'w', encoding='UTF8', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(fieldnames)
                for i in range(len(X_fuel_a)):
                    liste = [X_fuel_a[i], X_fuel_b[i], X_fuel_c[i], value[i]]
                    writer.writerow(liste)

    def plot_ternary_graph(self, a, b, c, value, pole_labels=['label1', 'label2', 'label3'], ncontours=8,
            showscale=True, colorscale='Portland', showmarkers=False, title=None):
        """Creates a ternary plot contour using plotly library

        Parameters
        ----------
        a : list
            list of proportions for the first axis variable
        b : list
            list of proportions for the second axis variable
        c : list
            list of proportions for the third axis variable
        value : list
            list of the variables for contour field corresponding to a, b, c
        pole_labels : list
            list of three labels to plot at the poles of the graph (Default value = ['label1', 'label2', 'label3'])
        ncontours : int
            number of division for the colorscale (Default value = 8)
        showscale :
            choose the color for the colorscale (Default value = True)
        colorscale : str
             (Default value = 'Portland')
        showmarkers : bool
            if True, show the points calculated, remark : points on boundaries do not appear (Default value = False)
        title : str
            if given, add a title to the figure (Default value = None)

        Warning
        -------
        Dependency with `plotly.figure_factory`

        """
        import plotly.figure_factory as ff

        return ff.create_ternary_contour(np.array([np.array(a), np.array(b), np.array(c)]), np.array(value),
                                         pole_labels=pole_labels, ncontours=ncontours, showscale=showscale,
                                         colorscale=colorscale, showmarkers=showmarkers, title=title)

    ###########################################################################
    # Features that are not yet functional and not priority for a 1.0 version #
    ###########################################################################

    def contour_plot(self, cases_list, mechanism, reactor_type, x, y, z, P=-1, T=-1, phi=-1,
                     show=True, save=False, compare_with=None, normalized=True, pcolor=False):
        """Plots a contour plot for 1 fixed parameters between P, T and phi

        Parameters
        ----------
        cases_list : list
            list of class :func:`~ARCANE.cases.Case` object
        mechanism : class :func:`~ARCANE.mechanisms.Mechanism` object
            class :func:`~ARCANE.mechanisms.Mechanism` object
        reactor_type : str
            type of reactor, `ARCANE.kwdict.reactor_labels`, typically 'C0DV', 'C0DP', 'C1DP', 'C1DCFD', 'C1DF'
        x : `str`
            name of the variable wanted as an abscissa
        y : str
            name of the variable wanted as an ordinate
        z : str
            name of the variable wanted in colours
        P : float
            fixed pressure parameter (Default value = -1)
        T : float
            fixed temperature parameter (Default value = -1)
        phi : float
            fixed equivalence parameter (Default value = -1)
        show : bool
            if False, doesn't display the figure (Default value = True)
        save : bool
            if True, saves the figure in pdf format (Default value = False)
        compare_with :
            either a list of class :func:`~ARCANE.cases.Case` objects
            or a class :func:`~ARCANE.mechanisms.Mechanism` objects (Default value = None)
        normalized : bool
            if True, the error is normalized by the value of the reference mechanism (Default value = True)
        pcolor : matplotlib.pyplot.pcolor
            plots as a pcolor, meaning, values without interpolation (Default value = False)

        See Also
        --------
        matplotlib.pyplot.pcolor

        """

        if isinstance(mechanism, list):
            if len(mechanism) > 1:
                logger.warning('WARNING First mechanisms of the list will be taken')
            mechanism = mechanism[0]

        number_of_input = 0
        fixed_variable = ['P', 'T', 'phi']
        if P > -1:
            number_of_input += 1
            fixed_variable = 'P'
        if T > -1:
            number_of_input += 1
            fixed_variable = 'T'
        if phi > -1:
            number_of_input += 1
            fixed_variable = 'phi'

        # Consistency checks

        if number_of_input < 1:
            logger.error('ERROR: You need to fix only 1 parameter for a contour plot')
            logger.error('You can specify P, T or phi with those keywords')
            sys.exit()

        selected_cases = tools.get_case_by_state(cases_list, reactor_type=reactor_type, P=P, T=T, phi=phi)
        if len(selected_cases) == 0:
            logger.error('ERROR: No matching cases for parametric plot')
            sys.exit()

        abscissa = []
        ordinate = []
        colours_new = []

        # Extracting quantity and method from string
        splitted = z.split(' ')
        last_split = splitted[-1]

        quantity = z
        method = ''

        # Checking if last word is a method
        for name in kwdict.methods:
            if last_split in kwdict.methods[name]:
                method = last_split
                quantity = ' '.join(splitted[:-1])

        for case in selected_cases:

            if fixed_variable == 'P':
                if x in kwdict.names['T']:
                    abscissa.append(float(case.temperature))
                    ordinate.append(float(case.phi))
                else:
                    abscissa.append(float(case.phi))
                    ordinate.append(float(case.temperature))
            elif fixed_variable == 'T':
                if x in kwdict.names['P']:
                    abscissa.append(float(case.pressure))
                    ordinate.append(float(case.phi))
                else:
                    abscissa.append(float(case.phi))
                    ordinate.append(float(case.pressure))
            elif fixed_variable == 'phi':
                if x in kwdict.names['T']:
                    abscissa.append(float(case.phi))
                    ordinate.append(float(case.pressure))
                else:
                    abscissa.append(float(case.phi))
                    ordinate.append(float(case.pressure))

            colours_new.append(case.extract_quantity(z, mechanism=mechanism))

        # Definition for comparison

        if compare_with:

            colours_ref = []

            if type(compare_with) == type(cases_list):

                selected_cases_ref = tools.get_case_by_state(compare_with, reactor_type=reactor_type, P=P, T=T, phi=phi)
                reference_mechanism = mechanism

            elif type(compare_with) == type(mechanism):

                selected_cases_ref = selected_cases
                reference_mechanism = compare_with

            for case in selected_cases_ref:

                if fixed_variable == 'P':
                    if x in kwdict.names['T']:
                        abscissa.append(float(case.temperature))
                        ordinate.append(float(case.phi))
                    else:
                        abscissa.append(float(case.phi))
                        ordinate.append(float(case.temperature))
                elif fixed_variable == 'T':
                    if x in kwdict.names['P']:
                        abscissa.append(float(case.pressure))
                        ordinate.append(float(case.phi))
                    else:
                        abscissa.append(float(case.phi))
                        ordinate.append(float(case.pressure))
                elif fixed_variable == 'phi':
                    if x in kwdict.names['T']:
                        abscissa.append(float(case.phi))
                        ordinate.append(float(case.pressure))
                    else:
                        abscissa.append(float(case.phi))
                        ordinate.append(float(case.pressure))

                colours_ref.append(case.extract_quantity(z, mechanism=reference_mechanism))

            if normalized:
                colours = [(c_ref - c_new) / c_ref for c_ref, c_new in zip(colours_ref, colours_new)]
            else:
                colours = [c_ref / c_new for c_ref, c_new in zip(colours_ref, colours_new)]

        else:
            colours = colours_new

        if quantity in mechanism.ctmech.species_names:
            quantity_name = 'Y'
        else:
            quantity_name = quantity

        colours_matrix = np.zeros([len(set(abscissa)), len(set(ordinate))])

        absci = []
        ordi = []
        for index, colour in enumerate(colours):
            x_val = abscissa[index]
            y_val = ordinate[index]

            if x_val not in absci:
                absci.append(x_val)

            if y_val not in ordi:
                ordi.append(y_val)

            colours_matrix[absci.index(x_val), ordi.index(y_val)] = colour

        if x == '1000/T':
            absci = [1000 / T for T in absci]
        elif y == '1000/T':
            ordi = [1000 / T for T in ordi]

        ordinate, abscissa = np.meshgrid(ordi, absci)
        colours = colours_matrix

        max_val = np.max(colours)
        min_val = np.min(colours)
        if min_val < 0 and max_val > 0:
            min_bound = - max_val
            cmap = 'coolwarm'
        else:
            min_bound = min_val
            cmap = 'summer'

        if max_val < 1:
            fmt = '%1.3f'
        else:
            fmt = '%1.0f'

        parametric_variable_names = kwdict.names['P'] + kwdict.names['T'] + kwdict.names['phi'] + ['1000/T']

        # Checking if the keyword is correct
        if x not in parametric_variable_names:
            logger.error('The keyword for the abscissa is not a valid one: ' + x)
            logger.error('Here is a list of the ones you can use: ' + str(', ').join(parametric_variable_names))
            return

        if y not in parametric_variable_names:
            logger.error('The keyword for the abscissa is not a valid one:' + y)
            logger.error('Here is a list of the ones you can use: ' + str(', ').join(parametric_variable_names))
            return

        data_names = list(case.names_dictionary(mechanism).keys())
        data_names += kwdict.names['sl'] + kwdict.names['tig'] + kwdict.names['thickness']

        # Checking if the keyword is correct
        if quantity not in data_names:
            logger.error('The keyword for the ordinate quantity is not a valid one: ' + quantity)
            logger.error('Here is a list of the ones you can use: ' + ', '.join(data_names))
            return

        # Checking if the keyword is correct
        possible_methods = ['']
        for name in kwdict.methods:
            possible_methods += kwdict.methods[name]

        if method not in possible_methods:
            logger.error('The keyword for the ordinate method is not a valid one: ' + method)
            logger.error('Here is a list of the ones you can use: ' + ', '.join(possible_methods))
            return

        fig, ax = plt.subplots(figsize=(20, 10))

        if not pcolor:
            if z == 'Ignition delay time' and not compare_with:
                exponents = range(0, 7)[::-1]
                filled_levels = [i * 10 ** -exp for exp in exponents for i in range(1, 10)][:-9]
                levels = [i * 10 ** -exp for exp in exponents for i in [1, 5]][:-1]

                csf = ax.contourf(abscissa, ordinate, colours, levels=filled_levels, norm=colors.LogNorm(), cmap=cmap)
                cs = ax.contour(abscissa, ordinate, colours, levels=levels, norm=colors.LogNorm(), colors='k')
                ax.clabel(cs, inline_spacing=20, fontsize=10, fmt='%1.0e')
            else:
                csf = ax.contourf(abscissa, ordinate, colours, 50, cmap=cmap, vmin=min_bound)
                cs = ax.contour(abscissa, ordinate, colours, 10, colors='k', vmin=min_bound)
                ax.clabel(cs, inline_spacing=20, fontsize=10, fmt=fmt)
        else:
            if z == 'Ignition delay time' and not compare_with:
                exponents = range(0, 7)[::-1]
                filled_levels = [i * 10 ** -exp for exp in exponents for i in range(1, 10)][:-9]

                csf = ax.pcolor(abscissa, ordinate, colours, cmap=cmap, levels=filled_levels)
            else:
                csf = ax.pcolor(abscissa, ordinate, colours, cmap=cmap)

        cbar = fig.colorbar(csf)
        addon = ''
        if z.endswith('integral'):
            cbar.ax.set_ylabel(z + ' [' + kwdict.units['integral' + reactor_type] + ']')
        elif compare_with:
            cbar.ax.set_ylabel(z + ' difference')
            addon = ' difference'
        else:
            cbar.ax.set_ylabel(z + ' [' + kwdict.units[kwdict.get_full_name(quantity_name)] + ']')

        ax.set_xlabel(x + ' [' + kwdict.units[kwdict.get_full_name(x)] + ']')
        ax.set_ylabel(y + ' [' + kwdict.units[kwdict.get_full_name(y)] + ']')

        title = z + addon + ' for ' + kwdict.reactor_labels[reactor_type][0] + ' at '
        if fixed_variable[0] == 'P':
            title += 'P = ' + str(P / 1e5) + 'bar'
        if fixed_variable[0] == 'T':
            title += 'T = ' + str(T) + ' K'
        if fixed_variable[0] == 'phi':
            title += 'phi = ' + str(phi)

        ax.set_title(title)

        if save:
            path_to_save = self.graphdir + '/Contour' + addon[1:] \
                           + 'X' + x.replace(' ', '_').replace('/', '-') \
                           + 'Y' + y.replace(' ', '_').replace('/', '-') \
                           + 'Z' + z.replace(' ', '_')
            path_to_save += kwdict.reactor_labels[reactor_type][0].replace(' ', '_')
            if fixed_variable[0] == 'P':
                path_to_save += 'P' + str(P / 1e5)
            if fixed_variable[0] == 'T':
                path_to_save += 'T' + str(T)
            if fixed_variable[0] == 'phi':
                path_to_save += 'phi' + str(phi)
            path_to_save += '.pdf'
            plt.savefig(path_to_save)
            logger.debug('Figure saved to ' + path_to_save)

        if show:
            plt.show()

    def surface_plot(self, cases_list, mechanism, reactor_type, x, y, z, P=-1, T=-1, phi=-1, fuel='',
                     show=True, save=False, hold=False, legend=None):
        """Plots a parametric plot for 2 fixed parameters between P, T, phi and fuel

        Parameters
        ----------
        cases_list : list
            list of class :func:`~ARCANE.cases.Case` objects
        mechanism : list
            list of class :func:`~ARCANE.mechanisms.Mechanism` objects
        reactor_type : str
            type of reactor, `ARCANE.kwdict.reactor_labels`, typically 'C0DV', 'C0DP', 'C1DP', 'C1DCFD', 'C1DF'
        x : str
            name of the variable wanted as an abscissa
        y : str
            name of the variable wanted as second abscissa
        z : str
            name of the variable wanted as an ordinate
        P : float
            fixed pressure parameter (Default value = -1)
        T : float
            fixed temperature parameter (Default value = -1)
        phi : float
            fixed equivalence parameter (Default value = -1)
        fuel : dict
            fuel composition parameter (Default value = '')
        show : bool
            if False, doesn't display the figure (Default value = True)
        save : bool
            if True, saves the figure in pdf format (Default value = False)
        hold : bool
            keep the previous graphs in memory (Default value = False)
        legend : str
            gives a non default label to the plot

        """

        if isinstance(mechanism, list):
            mechanism = mechanism[0]
        else:
            mechanism = mechanism

        number_of_input = 0
        parametric_variable_names = ['Pressure', 'Temperature', '1000/T', 'Equivalence ratio',
                                     'Fuel composition', 'Fuel ratio']
        parametric_variable = ['P', 'T', 'phi', 'fuel']
        if P > -1:
            number_of_input += 1
            parametric_variable_names.remove('Pressure')
            parametric_variable.remove('P')
        if T > -1:
            number_of_input += 1
            parametric_variable_names.remove('Temperature')
            parametric_variable_names.remove('1000/T')
            parametric_variable.remove('T')
        if phi > -1:
            number_of_input += 1
            parametric_variable_names.remove('Equivalence ratio')
            parametric_variable.remove('phi')
        if fuel:
            number_of_input += 1
            parametric_variable_names.remove('Fuel composition')
            parametric_variable_names.remove('Fuel ratio')
            parametric_variable.remove('fuel')

        if x in 'Fuel composition':
            parametric_variable[1] = 'fuel compo'
        elif x in 'Fuel ratio':
            parametric_variable[1] = 'fuel ratio'

        # Consistency checks
        fuels = []
        for case in cases_list:
            if case.fuel not in fuels:
                fuels.append(case.fuel)

        if not fuel and number_of_input < 2 and len(fuels) > 1:
            logger.warning("The fuel composition needs to be fixed as several are available :")
            quit(', '.join(fuels))

        if number_of_input < 2:
            logger.warning('WARNING: You need to fix 2 parameters for a parametric plot')
            quit('You can specify P, T, phi or fuel with those keywords')

        if x not in parametric_variable_names:
            quit('WARNING: x does not correspond to the variable parameter (' + x + ' and ' + parametric_variable_names[
                0] + ')')

        selected_cases = tools.get_case_by_state(cases_list, reactor_type=reactor_type, P=P, T=T, phi=phi, fuel=fuel)
        if len(selected_cases) == 0:
            quit('WARNING: No matching cases for parametric plot')

        # Plotting loop
        # for index_mech, mechanism in enumerate(mechanisms_list):
        abscissa_x = []
        abscissa_y = []
        ordinate = []

        index = 0
        secondary_fuel = ''
        while not secondary_fuel:
            fuel = cases_list[index].myid.split('T', 1)[1]
            fuel = fuel[4:]
            fuel = fuel.split('%', 1)
            primary_fuel = fuel[0]
            secondary_fuel = fuel[1][3:].split('%', 1)[0]
            index += 1

        for case in selected_cases:

            case.axis = []
            data_s = self.extract_scalar_values(case, mechanism)

            if parametric_variable[0] == 'P':
                abscissa_x.append(float(case.pressure))
            elif parametric_variable[0] == 'T':
                abscissa_x.append(case.temperature)
            elif parametric_variable[0] == 'phi':
                abscissa_x.append(float(case.phi))
            elif parametric_variable[0] == 'fuel compo':
                if primary_fuel in case.fuel:
                    abscissa_x.append(case.fuel[primary_fuel])
                else:
                    abscissa_x.append(0)
            elif parametric_variable[0] == 'fuel ratio':
                if primary_fuel in case.fuel:
                    if secondary_fuel in case.fuel:
                        abscissa_x.append(case.fuel[primary_fuel] / case.fuel[secondary_fuel])
                    else:
                        abscissa_x.append(abscissa_x[-1] + 10)
                else:
                    abscissa_x.append(0)

            if parametric_variable[1] == 'P':
                abscissa_y.append(float(case.pressure))
            elif parametric_variable[1] == 'T':
                abscissa_y.append(case.temperature)
            elif parametric_variable[1] == 'phi':
                abscissa_y.append(float(case.phi))
            elif parametric_variable[1] == 'fuel compo':
                if primary_fuel in case.fuel:
                    abscissa_y.append(case.fuel[primary_fuel])
                else:
                    abscissa_y.append(0)
            elif parametric_variable[1] == 'fuel ratio':
                if primary_fuel in case.fuel:
                    if secondary_fuel in case.fuel:
                        abscissa_y.append(case.fuel[primary_fuel] / case.fuel[secondary_fuel])
                    else:
                        abscissa_y.append(abscissa_y[-1] + 10)
                else:
                    abscissa_y.append(0)

            ordinate.append(float(data_s[self.available_scalars_dict[z]]))

        ordinate_matrix = np.zeros([len(set(abscissa_y)), len(set(abscissa_x))])

        absci_x = []
        absci_y = []

        for index, value in enumerate(ordinate):

            x_val = abscissa_x[index]
            y_val = abscissa_y[index]

            if x_val not in absci_x:
                absci_x.append(x_val)

            if y_val not in absci_y:
                absci_y.append(y_val)

            ordinate_matrix[absci_y.index(y_val), absci_x.index(x_val)] = value

        if x == '1000/T':
            abscissa_x = [1000 / T for T in abscissa_x]

        parametric_variable_names = ['Pressure', 'Temperature', '1000/T', 'Equivalence ratio',
                                     'Fuel composition', 'Fuel ratio']
        if x not in parametric_variable_names \
                or y not in parametric_variable_names \
                or z not in self.available_scalars_variables:
            logger.warning('WARNING: Wrong keywords for x or y')
            logger.warning('Available keywords for x are: ' + str(', ').join(parametric_variable_names))
            logger.warning('Available keywords for y are: ' + str(', ').join(self.available_scalars_variables))
            quit('Given: x = ' + x + ' and ' + 'y = ' + y)

        if not legend:
            label = mechanism.nickname
        else:
            label = legend

        abscissa_x, abscissa_y = np.meshgrid(absci_x, absci_y)

        # Plotting
        # if not plt.fignum_exists(1):
        fig = plt.figure(figsize=(20, 10))
        ax = fig.gca(projection='3d')

        ax.plot_surface(abscissa_x, abscissa_y, ordinate_matrix, label=label, cmap=cm.coolwarm)

        if y == 'Ignition delay time':
            plt.yscale('log')

        if x in 'Fuel ratio':
            ax.set_ylabel(x + ' $X_{' + primary_fuel + '}$ / $X_{' + secondary_fuel + '}$ [' + kwdict.units[x] + ']')
        elif x in 'Fuel composition':
            ax.set_ylabel(x + ' $X_{' + primary_fuel + '}$ [' + kwdict.units[x] + ']')
        else:
            ax.set_ylabel(x + ' [' + kwdict.units[x] + ']')

        if y in 'Fuel ratio':
            ax.set_xlabel(y + ' $X_{' + primary_fuel + '}$ / $X_{' + secondary_fuel + '}$ [' + kwdict.units[y] + ']')
        elif y in 'Fuel composition':
            ax.set_xlabel(y + ' $X_{' + primary_fuel + '}$ [' + kwdict.units[y] + ']')
        else:
            ax.set_xlabel(y + ' [' + kwdict.units[y] + ']')

        if z.endswith('integral'):
            ax.set_zlabel(z + ' [' + kwdict.units['integral' + reactor_type] + ']')
        else:
            ax.set_zlabel(z + ' [' + kwdict.units[z] + ']')

        title = kwdict.reactor_labels[reactor_type][0] + ' at '
        if parametric_variable[0] == 'P':
            title += 'P = ' + str(P / 1e5) + ' bar'
        if parametric_variable[0] == 'T':
            title += 'T = ' + str(T) + ' K'
        if parametric_variable[0] == 'phi':
            title += 'phi = ' + str(phi)
        if parametric_variable[1] == 'P':
            title += ' and P = ' + str(P / 1e5) + ' bar'
        if parametric_variable[1] == 'T':
            title += 'and T = ' + str(T) + ' K'
        if parametric_variable[1] == 'phi':
            title += 'and phi = ' + str(phi)

        if 'fuel' in parametric_variable[0]:
            title += 'P = ' + str(P / 1e5) + ' bar, ' + 'T = ' + str(T) + ' K and phi = ' + str(phi)

        plt.title(title)

        if save:  # and not hold:

            path_to_save = self.graphdir + '/Surface' + 'X' \
                           + x.replace(' ', '_').replace('/', '-') \
                           + 'Y' + y.replace(' ', '_')

            path_to_save += kwdict.reactor_labels[reactor_type][0].replace(' ', '_')
            if parametric_variable[0] == 'P':
                path_to_save += 'T' + str(T) + 'phi' + str(phi)
            if parametric_variable[0] == 'T':
                path_to_save += 'P' + str(P / 1e5) + 'phi' + str(phi)
            if parametric_variable[0] == 'phi':
                path_to_save += 'P' + str(P / 1e5) + 'T' + str(T)
            path_to_save += '.pdf'
            plt.savefig(path_to_save, dpi=2000)
            logger.info('Figure saved to ' + path_to_save)

        if show and not hold:
            plt.rcParams["figure.figsize"] = [20, 10]
            plt.show()

        if not hold:
            plt.clf()
