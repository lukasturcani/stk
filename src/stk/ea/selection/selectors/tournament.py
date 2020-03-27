
class Tournament(_BatchingSelector, Selector):
    """
    Yields batches of molecules through tournament selection.

    In tournament selection, a random number of batches is chosen from
    the population undergo a competition. In each competition, the
    batch with the highest fitness value is yielded. This is repeated
    until `num_batches` are yielded.

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        tournament = stk.Tournament(
            num_batches=5,
            batch_size=1
        )

        # Select the molecules.
        for selected, in tournament.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)

    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
        random_seed=None,
    ):
        """
        Initialize a :class:`.Tournament` instance.

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

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._generator = np.random.RandomState(random_seed)
        if num_batches is None:
            num_batches = float('inf')

        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._batch_size = batch_size
        self._num_batches = num_batches
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        # The tournament can only take place if there is more than 1
        # batch.
        while (
            len(batches) > 1 and yielded.get_num() < self._num_batches
        ):
            tournament_size = self._generator.randint(
                low=2,
                high=len(batches)+1
            )
            competitors = self._generator.choice(
                a=batches,
                size=tournament_size,
                replace=False
            )
            yield max(competitors)

            if not self._duplicate_mols:
                batches = filter(yielded.has_no_yielded_mols, batches)
            if not self._duplicate_batches:
                batches = filter(yielded.is_unyielded_batch, batches)
            if not self._duplicate_mols or not self._duplicate_batches:
                batches = tuple(batches)


