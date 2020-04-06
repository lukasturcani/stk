"""
Roulette
========

"""

from .selector import Selector


class Roulette(Selector):
    """
    Uses roulette selection to select batches of molecules.

    In roulette selection the probability a batch is selected
    is given by its fitness. If the total fitness is the sum of all
    fitness values, the chance a batch is selected is given
    by::

        p = batch fitness / total fitness,

    where ``p`` is the probability of selection and the batch
    fitness is the sum of all fitness values of molecules in the
    batch [#]_.

    References
    ----------
    .. [#] http://tinyurl.com/csc3djm

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        roulette = stk.Roulette()

        # Select the molecules.
        for selected, in roulette.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)

    Yielding multiple molecules at once. For example, if molecules need
    to be selected for crossover

    .. code-block:: python

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        roulette = stk.Roulette(batch_size=2)

        # Select the molecules.
        for selected in roulette.select(pop):
            # selected is a tuple of length 2, holding the selected
            # molecules. You can do stuff with the selected molecules
            # Like apply crossover operations on them.
            offspring = list(crosser.cross(*selected))

    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
        random_seed=None
    ):
        """
        Initialize a :class:`Roulette` instance.

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

        if num_batches is None:
            num_batches = float('inf')

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._generator = np.random.RandomState(random_seed)
        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        while batches and yielded.get_num() < self._num_batches:
            total = sum(batch.get_fitness() for batch in batches)
            weights = [
                batch.get_fitness() / total for batch in batches
            ]
            yield self._generator.choice(batches, p=weights)

            if not self._duplicate_mols:
                batches = filter(yielded.has_no_yielded_mols, batches)
            if not self._duplicate_batches:
                batches = filter(yielded.is_unyielded_batch, batches)
            if not self._duplicate_mols or not self._duplicate_batches:
                batches = tuple(batches)


