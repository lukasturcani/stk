"""
Random Crosser
==============

"""


import numpy as np


class RandomCrosser:
    """
    Use some other crosser at random.

    Examples
    --------
    *Use One of Several Crossers at Random*

    .. code-block:: python

        import stk

        crosser = stk.RandomCrosser(
            crossers=(
                # Assume each crosser crosses in a different way.
                stk.GeneticRecombination(...),
                stk.GeneticRecombination(...),
            ),
        )
        # Use one of the component crossers at random.
        crossover_record1 = crosser.cross(record)
        # A different crosser may get selected at random the second,
        # third, etc, time.
        crossover_record2 = crosser.cross(record)

    """

    def __init__(self, crossers, weights=None, random_seed=None):
        """
        Initialize a :class:`.RandomMutator` instance.

        Parameters
        ----------
        crossers : :class:`tuple`
            Holds instances which have :meth:`.cross` method. The
            :meth:`.cross` method must return an instance of
            :class:`.CrossoverRecord`.

        weights : :class:`tuple` of :class:`float`, optional
            For each mutator, the probability that it will be chosen
            whenever :meth:`.mutate` is called.
            If ``None`` all `mutators` will have equal chance of being
            selected.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._crossers = crossers
        self._weights = weights
        self._generator = np.random.RandomState(random_seed)

    def cross(self, records):
        """
        Cross `records`.

        Parameters
        ----------
        records : :class:`iterable` of :class:`.MoleculeRecord`
            The molecule records on which a crossover operation is
            performed.

        Yields
        -------
        :class:`.CrossoverRecord`
            A record of a crossover operation.

        """

        crosser = self._generator.choice(
            a=self._crossers,
            p=self._weights,
        )
        return crosser.cross(records)
