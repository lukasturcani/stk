
class Worst(_BatchingSelector, Selector):
    """
    Selects batches of molecules, lowest fitness value first.

    Examples
    --------
    Select the worst 5 batches of size 3

    .. code-block:: python

        import stk

        population = stk.Population(...)
        worst = stk.Worst(num_batches=5, batch_size=3)
        for batch in worst.select(population):
            # Do stuff with batch.


    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
    ):
        """
        Initialize a :class:`.Worst` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.
            Duplicate batches can occur if the same molecule is found
            multiple times in a population.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        """

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        batches = sorted(batches)

        if not self._duplicate_mols:
            batches = filter(yielded.has_no_yielded_mols, batches)

        if not self._duplicate_batches:
            batches = filter(yielded.is_unyielded_batch, batches)

        yield from it.islice(batches, self._num_batches)


