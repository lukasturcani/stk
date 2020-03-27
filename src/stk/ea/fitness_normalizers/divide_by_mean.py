
class DivideByMean(_FilteringNormalizer, FitnessNormalizer):
    """
    Divides fitness values by the population mean.

    While this function can be used if the fitness value of each
    :class:`.Molecule` in the :class:`.EAPopulation` is a single
    number, it is most useful when the fitness values is a
    :class:`list` of numbers. In this case, it is necessary to somehow
    combine the numbers so that a single fitness value is produced.
    For example, take a fitness value which is the vector holding the
    properties ``[energy, diameter, num_atoms]``. For a given molecule
    these numbers may be something like ``[200,000, 12, 140]``. If we
    were to sum these numbers, the energy term would dominate the final
    fitness value. In order to combine these numbers we can divide them
    by the population averages. For example, if the average energy
    of molecules in the population is ``300,000`` the average diameter
    is ``10`` and the average number of atoms is ``70`` then the
    fitness vector would be scaled to ``[0.5, 1.2, 2]``. These
    numbers are now of a similar magnitude and can be summed to give a
    reasonable value. After division , each value represents how
    much better than the population average each property value is.
    In essence we have removed the units from each parameter.

    Examples
    --------
    Scale fitness values

    .. code-block:: python

        import stk

        mean_scaler = stk.DivideByMean()

        # Normalize the fitness values.
        # Assume the fitness values are
        # {mol1: 1, mol2: 2, mol3: 3}
        normalized = mean_scaler.normalize(pop)

        # normalized is
        # {mol1: 0.5, mol2: 1, mol3: 1.5}


    Scale fitness vectors

    .. code-block:: python

        # Create the normalizer.
        # mean_scaler = DivideByMean()

        # Normalize the fitness values.
        # Assume the fitness values are
        # {mol1: [1, 10, 100], mol2: [2, 20, 100], mol3: [3, 30, 100]}.
        normalized = mean_scaler.normalize(pop)

        # normalized is
        # {
        #     mol1: [0.5, 0.5, 0.5],
        #     mol2: [1, 1, 1],
        #     mol3: [1.5, 1.5, 1.5]
        # }.

    """

    def __init__(self, filter=lambda population, mol: True):
        """
        Initialize a :class:`.DivideByMean` instance.

        Parameters
        ----------
        filter : :class:`callable`, optional
            Takes a two parameters, first is a :class:`.EAPopulation`
            and the second is a :class:`.Molecule`, and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.
            The :class:`.EAPopulation` on which :meth:`normalize` is
            called is passed as the first argument while the second
            argument will be passed every :class:`.Molecule` in it.

        """

        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        # filtered gets iterated through multiple times.
        filtered = list(filtered)
        mean = np.mean(
            a=[fitness_values[mol] for mol in filtered],
            axis=0,
        )
        logger.debug(f'Means used in DivideByMean: {mean}')

        for mol in filtered:
            # Use divide here so that both lists and arrays work
            # properly.
            yield mol, np.divide(fitness_values[mol], mean)


