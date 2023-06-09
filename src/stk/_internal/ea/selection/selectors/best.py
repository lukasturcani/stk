"""
Best
====

"""

import itertools as it

from stk.molecular import Inchi

from .selector import Selector


class Best(Selector):
    """
    Selects batches of molecules, highest fitness value first.

    Examples
    --------
    *Yielding Single Molecule Batches*

    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation.

    .. testcode:: yielding-single-molecule-batches

        import stk

        # Make the selector.
        best = stk.Best()

        population = (
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
            ).with_fitness_value(1),
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
            ).with_fitness_value(2)
        )

        # Select the molecules.
        for selected, in best.select(population):
            # Do stuff with each selected molecule.
            pass

    *Yielding Batches Holding Multiple Molecules*

    Yielding multiple molecules at once. For example, if molecules need
    to be selected for crossover.

    .. testcode:: yielding-batches-holding-multiple-molecules

        import stk

        # Make the selector.
        best = stk.Best(batch_size=2)

        population = (
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
            ).with_fitness_value(1),
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
            ).with_fitness_value(2)
        )

        # Select the molecules.
        for selected1, selected2 in best.select(population):
            # Do stuff with the selected molecules.
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
    ):
        """
        Initialize a :class:`.Best` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecule records yielded at once.

        duplicate_molecules : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.
            Duplicate batches can occur if the same molecule is found
            multiple times in a population.

        key_maker : :class:`.MoleculeKeyMaker`, optional
            Used to get the keys of molecules. If two molecules have
            the same key, they are considered duplicates.

        fitness_modifier : :class:`callable`, optional
            Takes the `population` on which :meth:`.select` is called
            and returns a :class:`dict`, which maps records in the
            `population` to the fitness values the :class:`.Selector`
            should use. If ``None``, the regular fitness values of the
            records are used.

        """

        if fitness_modifier is None:
            fitness_modifier = self._get_fitness_values

        self._duplicate_molecules = duplicate_molecules
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        super().__init__(
            key_maker=key_maker,
            fitness_modifier=fitness_modifier,
        )

    def _select_from_batches(self, batches, yielded_batches):
        batches = sorted(batches, reverse=True)

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

        yield from it.islice(batches, self._num_batches)
