"""
A module for defining plotting functions.

"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from . import fitness_calculators

plt.switch_backend('agg')


def fitness_epp(pop,
                plot_name=False,
                dump_name=None,
                xlabel='Generation'):
    """
    Plots the min, max and avg fitness values of each subpopulation.

    Also saves the plot data in a csv file.

    Parameters
    ----------
    pop : :class:`.GAPopulation`
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
        ``.csv`` extension.

    xlabel : :class:`str`, optional
        The label on the x-axis.

    Returns
    -------
    None : :class:`NoneType`

    """

    sns.set(style='darkgrid')
    df = pd.DataFrame()
    for i, subpop in enumerate(pop.populations, 1):
        # Ignore failed molecules from EPP.
        clean_pop = [x.fitness for x in subpop if x.fitness != 0.0001]
        data = [
            {'Generation': i,
             'Fitness': max(clean_pop) if clean_pop else 1e-4,
             'Type': 'Max'},
            {'Generation': i,
             'Fitness': min(clean_pop) if clean_pop else 1e-4,
             'Type': 'Min'},
            {'Generation': i,
             'Fitness': np.mean(clean_pop) if clean_pop else 1e-4,
             'Type': 'Mean'}
        ]

        df = df.append(data, ignore_index=True)

    # Save the plot data.
    if plot_name and dump_name is None:
        basename, ext = os.path.splitext(plot_name)
        dump_name = basename + '.csv'
    df.to_csv(dump_name)

    if plot_name:
        fig = plt.figure(figsize=[8, 4.5])

        palette = sns.color_palette('deep')
        sns.scatterplot(
            x='Generation',
            y='Fitness',
            hue='Type',
            palette={'Max': palette[3],
                     'Min': palette[0],
                     'Mean': palette[2]},
            data=df)

        plt.legend(bbox_to_anchor=(1.15, 1), prop={'size': 9})
        plt.tight_layout()
        fig.savefig(plot_name, dpi=500)
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

    Also saves the plot data in a csv file.

    Parameters
    ----------
    pop : :class:`.GAPopulation`
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
        ``.csv`` extension.

    xlabel : :class:`str`, optional
        The label on the x-axis.

    Returns
    -------
    None : :class:`NoneType`

    """

    sns.set(style='darkgrid')
    func_name = pop.ga_tools.fitness.name
    fitness_func = vars(fitness_calculators)[func_name]

    # Exit if a function which does not have param_labels was used.
    if not hasattr(fitness_func, 'param_labels'):
        return

    df = pd.DataFrame()
    for gen, sp in enumerate(pop.populations):
        for i, param_label in enumerate(fitness_func.param_labels):
            param_vals = [m.progress_params[func_name][i] for m in sp]

            data = [
                {'Generation': gen,
                 param_label: max(param_vals) if sp else None,
                 'Type': 'Max'},
                {'Generation': gen,
                 param_label: min(param_vals) if sp else None,
                 'Type': 'Min'},
                {'Generation': gen,
                 param_label: np.mean(param_vals) if sp else None,
                 'Type': 'Mean'}]
            df = df.append(data, ignore_index=True)

    # Save the plot data.
    if plot_name and dump_name is None:
        basename, ext = os.path.splitext(plot_name)
        dump_name = basename + '.csv'
    df.to_csv(dump_name)

    palette = sns.color_palette('deep')
    for i, param_label in enumerate(fitness_func.param_labels):

        if plot_name:
            fig = plt.figure(figsize=[8, 4.5])

            sns.scatterplot(
                x='Generation',
                y=param_label,
                hue='Type',
                palette={'Max': palette[3],
                         'Min': palette[0],
                         'Mean': palette[2]},
                data=df)

            plt.legend(bbox_to_anchor=(1.15, 1), prop={'size': 9})
            plt.tight_layout()
            new_plot_name = str(i).join(os.path.splitext(plot_name))
            fig.savefig(new_plot_name, dpi=500)
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

    sns.set(style='darkgrid')
    fig = plt.figure()

    df = pd.DataFrame()
    for ind, selection_count in counter.items():
        label = f'{ind.name} - {ind.fitness}'
        data = {
            'Molecule: name - fitness value': label,
            'Number of times selected': selection_count,
            'Fitness': ind.fitness
        }
        df = df.append(data, ignore_index=True)

    df = df.sort_values(['Number of times selected', 'Fitness'],
                        ascending=[False, False])
    norm = plt.Normalize(df['Fitness'].min(), df['Fitness'].max())
    sm = plt.cm.ScalarMappable(cmap='magma_r', norm=norm)
    sm.set_array([])

    ax = sns.barplot(
                x='Molecule: name - fitness value',
                y='Number of times selected',
                hue='Fitness',
                palette='magma_r',
                dodge=False,
                data=df)
    ax.get_legend().remove()
    ax.figure.colorbar(sm).set_label('Fitness')
    plt.xticks(rotation=90)
    plt.tight_layout()
    fig.savefig(plot_name, dpi=fig.dpi)
    plt.close('all')
