"""
Random Mutator
==============

"""


import numpy as np


class RandomMutator:
    """
    Use some other mutator at random.

    Examples
    --------
    *Use One of Several Mutators at Random*

    .. code-block:: python

        import stk

        mutator = stk.RandomMutator(
            mutators=(
                stk.RandomBuildingBlock(...),
                stk.SimilarBuildingBlock(...),
                stk.RandomBuildingBlock(...),
            ),
        )
        # Use one of the component mutators at random.
        mutation_record = mutator.mutate(record)
        # A different mutator may get selected at random the second,
        # third, etc, time.
        mutation_record2 = mutator.mutate(record)

    """

    def __init__(self, mutators, weights=None, random_seed=None):
        """
        Initialize a :class:`.RandomMutator` instance.

        Parameters
        ----------
        mutators : :class:`tuple`
            Holds instances which have :meth:`.mutate` method. The
            :meth:`.mutate` method must return an instance of
            :class:`.MutationRecord`.

        weights : :class:`tuple` of :class:`float`, optional
            For each mutator, the probability that it will be chosen
            whenever :meth:`.mutate` is called.
            If ``None`` all `mutators` will have equal chance of being
            selected.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._mutators = mutators
        self._weights = weights
        self._generator = np.random.RandomState(random_seed)

    def mutate(self, record):
        """
        Return a mutant of `record`.

        Parameters
        ----------
        record : :class:`.MoleculeRecord`
            The molecule to be mutated.

        Returns
        -------
        :class:`.MutationRecord`
            A record of the mutation. The exact subclass of
            :class:`.MutationRecord` depends on which mutator was
            used.

        """

        mutator = self._generator.choice(
            a=self._mutators,
            p=self._weights,
        )
        return mutator.mutate(record)
