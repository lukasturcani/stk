"""
Above Average
=============

"""

import numpy as np
import itertools as it

from stk.molecular import Inchi
from .selector import Selector
from .utilities import Batcher


class AboveAverage(Selector):
    """
    Yields above average batches of molecules.

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        above_avg = stk.AboveAverage()

        # Select the molecules.
        for selected, in above_avg.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)

    Yielding multiple molecules at once. For example, if molecules need
    to be selected for crossover.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        above_avg = stk.AboveAverage(batch_size=2)

        # Select the molecules.
        for selected in above_avg.select(pop):
            # selected is a tuple of length 2, holding the selected
            # molecules. You can do stuff with the selected molecules
            # Like apply crossover operations on them.
            offspring = list(crosser.cross(*selected))

    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_records=True,
        duplicate_batches=True,
        key_maker=Inchi(),
        fitness_modifier=None,
    ):
        """
        Initialize an :class:`.AboveAverage` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecule records in each yielded
            :class:`.Batch`.

        duplicate_records : :class:`bool`, optional
            If ``True`` the same molecule record can be yielded in more
            than one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        key_maker : :class:`.MoleculeKeyMaker`, optional
            Used to get the keys of molecules, which are used to
            determine if two molecules are duplicates of each other.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        """

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._duplicate_records = duplicate_records
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        self._key_maker = key_maker
        self._fitness_modifier = fitness_modifier
        self._batcher = Batcher()

    def select(self, population):
        mean = np.mean([batch.get_fitness() for batch in batches])
        # Yield highest fitness batches first.
        batches = sorted(batches, reverse=True)
        # Yield only batches with a fitness larger than the mean.
        batches = it.takewhile(
            lambda batch: batch.get_fitness() > mean,
            batches
        )
        # Yield batches which are multiple times better than the mean
        # multiple times.
        batches = (
            batch
            for batch in batches
            for i in range(self._get_num_duplicates(batch, mean))
        )
        # If duplicate molecules are not allowed, allow only
        # batches with no yielded molecules.
        if not self._duplicate_mols:
            batches = filter(yielded.has_no_yielded_mols, batches)
        # If duplicate batches are not allowed, allow only
        # unyielded batches.
        if not self._duplicate_batches:
            batches = filter(yielded.is_unyielded_batch, batches)
        # Limit the number of yielded batches to _num_batches.
        yield from it.islice(batches, self._num_batches)

    def _get_num_duplicates(self, batch, mean):
        if self._duplicate_batches and self._duplicate_mols:
            return int(batch.get_fitness() // mean)
        return 1
