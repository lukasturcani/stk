
class ProgressPlotter(Plotter):
    """
    Plots how a property changes during a GA run.

    The produced plot will show the GA generations on the x axis and
    the min, mean and max values of an attribute on the y axis.

    Examples
    --------
    Plot how fitness value changes with GA generations

    .. code-block:: python

        import stk

        # progress has subpopulations where each subpopulation is a
        # generation of a GA.
        progress = stk.Population(...)

        # Make the plotter which plots the fitness change across
        # generations.
        plotter = stk.ProgressPlotter(
            filename='fitness_plot',
            property_fn=lambda mol: mol.fitness,
            y_label='Fitness'
        )
        plotter.plot(progress)

    Plot how the number of atoms changes with GA generations

    .. code-block:: python

        plotter = ProgressPlotter(
            filename='atom_plot',
            property_fn=lambda mol: len(mol.atoms),
            y_label='Number of Atoms'
        )
        plotter.plot(progress)

    """

    def __init__(
        self,
        filename,
        property_fn,
        y_label,
        progress_fn=None,
        filter=lambda progress, mol: True,
    ):
        """
        Initialize a :class:`ProgressPlotter` instance.

        Parameters
        ----------
        filename : :class:`str`
            The basename of the files. This means it should not include
            file extensions.

        property_fn : :class:`callable`
            A :class:`callable` which takes a :class:`.EAPopulation`
            and a :class:`.Molecule` object and returns a property
            value of that molecule, which is used for the plot.
            The :class:`callable` must return a valid value for each
            :class:`.Molecule` in the population.

        y_label : :class:`str`
            The y label for the produced graph.

        progress_fn : :class:`callable`, optional
            Takes the population passed to :meth:`plot` and excutes a
            computation on it. This may be useful if you want to
            apply a normalization to the fitness values in the
            progress population, for example.

        filter : :class:`callable`, optional
            Takes an :class:`.EAPopulation` and a :class:`.Molecule` as
            input and returns ``True`` or ``False``. Only molecules
            which return ``True`` will be plotted. Default is for all
            molecules to be plotted.

        """

        self._filename = filename
        self._property_fn = property_fn
        self._y_label = y_label
        self._filter = filter
        self._progress_fn = progress_fn

    def plot(self, progress):
        """
        Plot a progress plot.

        Parameters
        ----------
        progress : :class:`.Population`
            A :class:`.Population` where each generation of the GA is
            a subpopulation.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self._progress_fn is not None:
            self._progress_fn(progress)

        def filter_fn(mol):
            return self._filter(progress, mol)

        def property_fn(mol):
            return self._property_fn(progress, mol)

        sns.set(style='darkgrid')
        df = pd.DataFrame()
        for i, subpop in enumerate(progress.subpopulations, 1):
            filtered = filter(filter_fn, subpop)
            subpop_vals = list(map(property_fn, filtered))

            # If there are no values after filtering, don't plot
            # anything for the generation.
            if not subpop_vals:
                continue

            data = [
                {
                    'Generation': i,
                    self._y_label: max(subpop_vals),
                    'Type': 'Max'
                },
                {
                    'Generation': i,
                    self._y_label: min(subpop_vals),
                    'Type': 'Min'
                },
                {
                    'Generation': i,
                    self._y_label: np.mean(subpop_vals),
                    'Type': 'Mean'
                }
            ]

            df = df.append(data, ignore_index=True)

        # Save the plot data.
        df.to_csv(f'{self._filename}.csv')

        fig = plt.figure(figsize=[8, 4.5])

        palette = sns.color_palette('deep')
        # It's possible that all values were filtered out, and trying
        # to plot an empty dataframe would raise an exception.
        if len(df) != 0:
            sns.scatterplot(
                x='Generation',
                y=self._y_label,
                hue='Type',
                palette={
                    'Max': palette[3],
                    'Min': palette[0],
                    'Mean': palette[2]
                },
                data=df
            )
        # Set the length of the axes to account for all generations,
        # as its possible the first or last ones were not included
        # due to being filtered out.
        plt.xlim(0, len(progress.subpopulations)+1)

        plt.legend(bbox_to_anchor=(1.15, 1), prop={'size': 9})
        plt.tight_layout()
        fig.savefig(f'{self._filename}.png', dpi=500)
        plt.close('all')


