"""
Filter Molecules
================

"""

from .selector import Selector


class FilterMolecules(Selector):
    """
    Allows a :class:`.Selector` to select only some molecules.

    Examples
    --------
    *Using a Selection Algorithm on a Subset of Batches*

    .. testcode:: using-a-selection-algorithm-on-a-subset-of-batches

        import stk

        selector = stk.FilterMolecules(
            # Use only molecules in the top 5 batches of size 3.
            filter=stk.Best(num_batches=5, batch_size=3),
            # Select from those molecules using roulette selection.
            selector=stk.Roulette(num_batches=20, batch_size=3),
        )

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
            for i in range(10)
        )

        for batch in selector.select(population):
            # Do stuff with batch. All the molecules in the batch
            # belong to the top 5 batches of size 3. The batch
            # was selected using roulette selection.
            pass

    """

    def __init__(self, filter, selector):
        """
        Initialize a :class:`.FilterMolecules` instance.

        Parameters
        ----------
        filter : :class:`.Selector`
            Selects molecules which can be yielded by `selector`.

        selector : :class:`.Selector`
            Selects batches of molecules. The batches can only
            contain molecules yielded by `filter`.

        """

        self._filter = filter
        self._selector = selector

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None,
    ):
        valid_batches = self._filter.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
        valid_population = tuple(
            record
            for batch in valid_batches
            for record in batch
        )
        yield from self._selector.select(
            population=valid_population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
