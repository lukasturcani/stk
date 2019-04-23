"""

"""

from functools import wraps
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from collections import Counter

plt.switch_backend('agg')


class Plotter:
    """

    """

    def plot(self, progress):
        raise NotImplementedError()


class ProgressPlotter(Plotter):
    """

    """

    def __init__(self, attr, y_label, default):
        """

        """

        self.attr = attr
        self.y_label = y_label
        self.default = default

    def plot(self, progress):
        """

        """

        sns.set(style='darkgrid')
        df = pd.DataFrame()
        for i, subpop in enumerate(progress.populations, 1):
            # Ignore failed molecules from progress plot.
            pop = [
                getattr(x, self.attr)
                for x in subpop if getattr(x, self.attr) is not None
            ]
            data = [
                {'Generation': i,
                 self.y_label: max(pop) if pop else self.default,
                 'Type': 'Max'},
                {'Generation': i,
                 self.y_label: min(pop) if pop else self.default,
                 'Type': 'Min'},
                {'Generation': i,
                 self.y_label: np.mean(pop) if pop else self.default,
                 'Type': 'Mean'}
            ]

            df = df.append(data, ignore_index=True)

        # Save the plot data.
        df.to_csv(f'{self.attr}_progress.csv')

        fig = plt.figure(figsize=[8, 4.5])

        palette = sns.color_palette('deep')
        sns.scatterplot(
            x='Generation',
            y=self.y_label,
            hue='Type',
            palette={'Max': palette[3],
                     'Min': palette[0],
                     'Mean': palette[2]},
            data=df
        )

        plt.legend(bbox_to_anchor=(1.15, 1), prop={'size': 9})
        plt.tight_layout()
        fig.savefig(f'{self.attr}_progress.png', dpi=500)
        plt.close('all')


class SelectionPlotter(Plotter):
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

    def __init__(self, filename, selector):
        self.plots = 0
        self.filename = filename
        selector.select = self.update_counter(selector.select)

    def update_counter(self, select):
        """

        """

        @wraps(select)
        def inner(*args, **kwargs):
            selected = select(*args, **kwargs)
            counter = Counter()

            while True:
                try:
                    yielded = next(selected)
                except StopIteration:
                    self._plot(counter)
                    raise

                counter.update(yielded)
                yield yielded

        return inner

    def _plot(self, counter):
        self.plots += 1
        sns.set(style='darkgrid')
        fig = plt.figure()

        df = pd.DataFrame()
        for mol, selection_count in counter.items():
            label = f'{mol.name} - {mol.fitness}'
            data = {
                'Molecule: name - fitness value': label,
                'Number of times selected': selection_count,
                'Fitness': mol.fitness
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
        fig.savefig(f'{self.filename}_{self.plots}.png', dpi=fig.dpi)
        plt.close('all')

    def plot(self, progress):
        """

        """

        return
