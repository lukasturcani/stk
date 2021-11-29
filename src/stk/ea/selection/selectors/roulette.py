"""
Roulette
========

"""

import numpy as np

from stk.molecular import Inchi

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
    *Yielding Single Molecule Batches*

    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation

    .. testcode:: yielding-single-molecule-batches

        import stk

        # Make the selector.
        roulette = stk.Roulette(num_batches=5)

        population = tuple(
            stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ).with_fitness_value(i)
            for i in range(100)
        )

        # Select the molecules.
        for selected, in roulette.select(population):
            # Do stuff with each selected molecule.
            pass


    *Yielding Batches Holding Multiple Molecules*

    Yielding multiple molecules at once. For example, if molecules need
    to be selected for crossover

    .. testcode:: yielding-batches-holding-multiple-molecules

        import stk

        # Make the selector.
        roulette = stk.Roulette(num_batches=5, batch_size=2)

        population = tuple(
            stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ).with_fitness_value(i)
            for i in range(100)
        )

        # Select the molecules.
        for selected1, selected2 in roulette.select(population):
            # Do stuff to the molecules.
            pass

    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_molecules=True,
        duplicate_batches=True,
        key_maker=Inchi(),
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

        duplicate_molecules : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        key_maker : :class:`.MoleculeKeyMaker`, optional
            Used to get the keys of molecules. If two molecules have
            the same key, they are considered duplicates.

        fitness_modifier : :class:`callable`, optional
            Takes the `population` on which :meth:`.select` is called
            and returns a :class:`dict`, which maps records in the
            `population` to the fitness values the :class:`.Selector`
            should use. If ``None``, the regular fitness values of the
            records are used.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        if num_batches is None:
            num_batches = float('inf')

        if fitness_modifier is None:
            fitness_modifier = self._get_fitness_values

        self._generator = np.random.RandomState(random_seed)
        self._duplicate_molecules = duplicate_molecules
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        super().__init__(
            key_maker=key_maker,
            fitness_modifier=fitness_modifier,
        )

    def _select_from_batches(self, batches, yielded_batches):
        while (
            batches and yielded_batches.get_num() < self._num_batches
        ):
            total = sum(batch.get_fitness_value() for batch in batches)
            weights = [
                batch.get_fitness_value() / total for batch in batches
            ]
            yield self._generator.choice(batches, p=weights)

            if not self._duplicate_molecules:
                batches = filter(
                    yielded_batches.has_no_yielded_molecules,
                    batches,
                )
            if not self._duplicate_batches:
                batches = filter(
                    yielded_batches.is_unyielded_batch,
                    batches,
                )
            if (
                not self._duplicate_molecules
                or not self._duplicate_batches
            ):
                batches = tuple(batches)
