"""
A module for defining plotting functions.

"""

import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from operator import attrgetter

from .fitness import *

def epp(progress, plot_name, xlabel='Generation'):
    """
    Plots all the EPPs.

    Parameters
    ----------
    progress : Population
        A population which holds all of the generations in a separate
        subpopulation.

    plot_name : str
        The full path of the .png file to which the plots should be
        saved. Values may be attached to the end of this name if
        multiple plots are plotted.

    xlabel : str (default = 'Generation')
        The label on the x-axis of the EPP plot.

    Modifies
    --------
    Creates plots at the path of `plot_name`.

    If `norm_func` is provided the data in `progress` is replaced with
    normalized values.

    Returns
    -------
    None : NoneType

    """

    # ``True`` if the fitness function defined `progress_params`. In
    # this case plot the EPP for each of the `progress_params`.
    if any(hasattr(ind, 'progress_params') for ind in progress):
        parameter_epp(progress, plot_name)

    # Renormalize all the fitness values across the entire GA run if
    # a Normalization object was provided. If previously the values
    # progress.maxs (and mins etc) were lists, they will be converted
    # to ints/floats.
    if norm:
        progress.normalize(norm)

    # Plot the EPP of the fitness values.
    fitness_epp(progress, plot_name, xlabel)

def fitness_epp(pop, plot_name, xlabel='Generation'):
    """
    Plots the min, max and avg fitness values of each subpopulation.

    Parameters
    ----------
    pop : Population
        A population containing any number of subpopulations and no
        direct members.

    plot_name : str
        The full path of where the plot should be saved.

    xlabel : str (default = 'Generation')
        The label on the x-axis of the EPP plot.

    Returns
    -------
    None : NoneType

    """

    xvals = []
    maxs = []
    means = []
    mins = []

    for i, subpop in enumerate(pop.populations, 1):
        xvals.append(i)
        maxs.append(max(x.fitness for x in subpop))
        means.append(subpop.mean(lambda x : x.fitness))
        mins.append(min(x.fitness for x in subpop))

    fig = plt.figure()
    plt.xlabel(xlabel)
    plt.ylabel('Fitness')
    plt.scatter(xvals, maxs, color='red', marker='x', label='max')
    plt.scatter(xvals, means, color='green', marker='x', label='mean')
    plt.scatter(xvals, mins, color='blue', marker='x', label='min')
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.savefig(plot_name, dpi=1000,
                bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close('all')

def parameter_epp(pop, plot_name):
    """
    Plots the progress_params values across subpopulations.

    For each progress_param a plot will be produced. Each will show
    subpopulations on the x-axis and the min, max and avg values
    of that progress param on the y axis.

    Parameters
    ----------
    pop : Population
        A population holding any number of subpopulations and no
        direct members.

    plot_name : str
        The full path of where the plots should be saved.

    Returns
    -------
    None : NoneType

    """

    func_name = pop.ga_tools.input.fitness_func.name
    fitness_func = globals()[func_name]

    min_params = []
    max_params = []
    mean_params = []
    xvals = list(range(1, len(pop.populations)+1 ))

    for sp in pop.populations:

        unscaled_var_mat = np.matrix([x.progress_params for x in sp])

        max_params.append(np.max(unscaled_var_mat,
                                    axis=0).tolist()[0])
        min_params.append(np.min(unscaled_var_mat,
                                    axis=0).tolist()[0])
        mean_params.append(np.mean(unscaled_var_mat,
                                    axis=0).tolist()[0])

    for x in range(len(min_params[0])):
        fig = plt.figure()
        plt.xlabel('Population')
        plt.ylabel('Unscaled ' + fitness_func.param_labels[x])
        plt.title('Population Comparison', fontsize=18)

        y_mean = [array[x] for array in mean_params]
        y_max = [array[x] for array in max_params]
        y_min = [array[x] for array in min_params]

        plt.scatter(xvals, y_mean,
                    color='green', marker='x', label='mean')
        plt.scatter(xvals, y_min,
                    color='blue', marker='x', label='min')
        plt.scatter(xvals, y_max,
                    color='red', marker='x', label='max')
        lgd=plt.legend(bbox_to_anchor=(1.05, 1), loc=2,
                        borderaxespad=0.)

        new_plot_name = str(x).join(os.path.splitext(plot_name))

        fig.savefig(new_plot_name, dpi=1000,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close('all')

def plot_counter(counter, plot_name):
    """
    Saves a .png file holding a plot of `counter`.

    The counter should hold the number of times a certain population
    member was selected.

    Parameters
    ----------
    counter : Counter
        A counter of which members of a population were picked.

    plot_name : str
        The full path of the .png where the plot is to be saved.

    Returns
    -------
    None : NoneType

    """

    fig = plt.figure()
    x_vals = list(range(1, len(counter.items())+1))
    y_vals = []
    labels = []

    for ind, value in sorted(counter.items(), reverse=True):
        y_vals.append(value)
        labels.append(ind.fitness)

    plt.bar(x_vals, y_vals, color='blue')
    plt.xticks([x+0.5 for x in x_vals], labels, rotation='vertical')
    plt.tight_layout()
    fig.savefig(plot_name, dpi=fig.dpi)
    plt.close('all')
