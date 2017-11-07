"""
A module for defining plotting functions.

"""

import os
import matplotlib.pyplot as plt
import numpy as np

from . import fitness

plt.switch_backend('agg')


def fitness_epp(pop,
                plot_name=False,
                dump_name=None,
                xlabel='Generation'):
    """
    Plots the min, max and avg fitness values of each subpopulation.

    Also saves the plot data. The min, mean and max values for each
    subpopulation of `pop` are placed in an array:

        >>> dmp_array = np.array([mins, means, maxs])

    which is then dumped to the file `dump_name`.

    The members with a fitness of ``0.0001`` (indicating a failed
    fitness calculation) are excluded from this analysis.

    Parameters
    ----------
    pop : :class:`.Population`
        A population containing any number of subpopulations and no
        direct members. Each subpopulation represents a generation of
        the GA.

    plot_name : :class:`str` or :class:`bool`, optional
        The full path to where the plot should be saved. If ``False``
        then no plot is made. Useful if you only want the dump file and
        not graph.

    dump_name : :class:`str`, optional
        The path to the file where the data to make the graph is
        dumped. If ``None`` then `plot_name` is used but with a
        ``.dmp`` extension.

    xlabel : :class:`str`, optional
        The label on the x-axis.

    Returns
    -------
    None : :class:`NoneType`

    """

    xvals = []
    maxs = []
    means = []
    mins = []

    for i, subpop in enumerate(pop.populations, 1):
        xvals.append(i)
        if len(subpop) == 0:
            maxs.append(0)
            means.append(0)
            mins.append(0)
        else:
            maxs.append(max(x.fitness for x in subpop))
            means.append(np.mean([x.fitness for x in subpop]))
            mins.append(min(x.fitness for x in subpop))

    # Save the plot data.
    if plot_name and dump_name is None:
        basename, ext = os.path.splitext(plot_name)
        dump_name = basename + '.dmp'

    np.array([mins, means, maxs]).dump(dump_name)

    if plot_name:
        fig = plt.figure()
        axes = plt.gca()
        plt.xlabel(xlabel)
        plt.ylabel('Fitness')
        plt.scatter(xvals, maxs, color='red', marker='o', label='max',
                    alpha=0.5)
        plt.scatter(xvals, means, color='green', marker='o',
                    label='mean', alpha=0.5)
        plt.scatter(xvals, mins, color='blue', marker='o', label='min',
                    alpha=0.5)
        frame = (min(mins) + max(maxs))/50
        axes.set_ylim([min(mins) - frame, max(maxs) + frame])

        lgd = plt.legend(bbox_to_anchor=(1.05, 1),
                         loc=2, borderaxespad=0.)
        fig.savefig(plot_name, dpi=1000,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close('all')


def parameter_epp(pop,
                  plot_name=False,
                  dump_name=None,
                  xlabel='Generation'):
    """
    Plots the :attr:`.MacroMolecule.progress_params` values.

    For each element in :attr:`~.MacroMolecule.progress_params` a
    separate graph will be produced. Each will show subpopulations on
    the x-axis and the min, max and avg values of that
    :attr:`~.MacroMolecule.progress_params` element on the y-axis.

    Also saves the plot data. The min, mean and max values for each
    subpopulation are placed in an array:

        >>> dmp_array = np.array([mins, means, maxs])

    which is then dumped to the file `dump_name`.

    Parameters
    ----------
    pop : :class:`.Population`
        A population containing any number of subpopulations and no
        direct members. Each subpopulation represents a generation of
        the GA.

    plot_name : :class:`str` or :class:`bool`, optional
        The full path to where the plot should be saved. If ``False``
        then no plot is made. Useful if you only want the dump file and
        not graph.

    dump_name : :class:`str`, optional
        The path to the file where the data to make the graph is
        dumped. If ``None`` then `plot_name` is used but with a
        ``.dmp`` extension.

    xlabel : :class:`str`, optional
        The label on the x-axis.

    Returns
    -------
    None : :class:`NoneType`

    """

    func_name = pop.ga_tools.fitness.name
    fitness_func = vars(fitness)[func_name]

    # Exit if a function which does not have param_labels was used.
    if not hasattr(fitness_func, 'param_labels'):
        return

    min_params = []
    max_params = []
    mean_params = []
    nparams = len(fitness_func.param_labels)

    for sp in pop.populations:

        if len(sp) == 0:
            min_params.append([None for x in range(nparams)])
            mean_params.append([None for x in range(nparams)])
            max_params.append([None for x in range(nparams)])
            continue

        p_mat = np.array([x.progress_params[func_name] for x in sp])

        # Each element of this list holds an array of all the valid
        # values of a particular progress_param.
        p_arrays = [p_mat[:, x] for x in range(nparams)]
        p_arrays = [[x for x in a if x is not None] for a in p_arrays]

        min_params.append([min(x) if len(x) > 0 else None
                           for x in p_arrays])

        mean_params.append([np.mean(x) if len(x) > 0 else None
                            for x in p_arrays])

        max_params.append([max(x) if len(x) > 0 else None
                           for x in p_arrays])

    for x in range(len(min_params[0])):

        y_mean = [array[x] for array in mean_params]
        y_max = [array[x] for array in max_params]
        y_min = [array[x] for array in min_params]

        # This loop checks if any y values are ``None`` in which case
        # they are removed and the corresponding x value as well.
        nones = set()
        for i, yval in enumerate(y_min):
            if yval is None:
                nones.add(i)
        xvals = list(range(1, len(pop.populations)+1))
        xvals = [xv for i, xv in enumerate(xvals) if i not in nones]
        y_mean = [yv for i, yv in enumerate(y_mean) if i not in nones]
        y_max = [yv for i, yv in enumerate(y_max) if i not in nones]
        y_min = [yv for i, yv in enumerate(y_min) if i not in nones]

        # Save the plot data.
        if plot_name and dump_name is None:
            basename, ext = os.path.splitext(plot_name)
            dump_name = basename+'.dmp'
        dname = dump_name.replace('.dmp', '{}.dmp'.format(x))
        np.array([y_min, y_mean, y_max]).dump(dname)

        if plot_name:
            fig = plt.figure()
            axes = plt.gca()
            plt.xlabel(xlabel)
            plt.ylabel('Unscaled ' + fitness_func.param_labels[x])
            plt.scatter(xvals, y_mean,
                        color='green', marker='o', label='mean')
            plt.scatter(xvals, y_min,
                        color='blue', marker='o', label='min')
            plt.scatter(xvals, y_max,
                        color='red', marker='o', label='max')
            frame = (min(y_min) + max(y_max))/50
            axes.set_ylim([min(y_min) - frame, max(y_max) + frame])
            lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2,
                             borderaxespad=0.)

            new_plot_name = str(x).join(os.path.splitext(plot_name))

            fig.savefig(new_plot_name, dpi=1000,
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close('all')


def plot_counter(counter, plot_name):
    """
    Saves a ``.png`` file holding a plot of `counter`.

    The counter should hold the number of times a certain population
    member was selected.

    Parameters
    ----------
    counter : :class:`collections.Counter`
        A counter of which members of a population were selected.

    plot_name : :class:`str`
        The full path of the ``.png`` where the plot is to be saved.

    Returns
    -------
    None : :class:`NoneType`

    """

    fig = plt.figure()
    x_vals = list(range(1, len(counter.items())+1))
    y_vals = []
    labels = []

    for ind, value in sorted(counter.items(), reverse=True):
        y_vals.append(value)
        labels.append(ind.name + ' - ' + str(ind.fitness))

    plt.bar(x_vals, y_vals, color='blue')
    plt.xlabel('Individuals, denoted by fitness value')
    plt.ylabel('Number of times selected')
    plt.xticks([x for x in x_vals], labels, rotation='vertical')
    plt.tight_layout()
    fig.savefig(plot_name, dpi=fig.dpi)
    plt.close('all')
