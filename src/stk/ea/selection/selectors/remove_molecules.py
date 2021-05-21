"""
Remove Molecules
================

"""

from .selector import Selector


class RemoveMolecules(Selector):
    """
    Prevents a :class:`.Selector` from selecting some molecules.

    Examples
    --------
    *Removing Molecules From Selection*

    .. testcode:: removing-molecules-from-selection

        import stk

        selector = stk.RemoveMolecules(
            # Do not select any molecules from the top 5 batches
            # of size 3.
            remover=stk.Best(num_batches=5, batch_size=3),
            # Select the ret of the molecules using roulette.
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
            # Do stuff with batch. The batch is guaranteed not to
            # contain any molecules which are found in the best 5
            # batches of size 3.
            pass

    """

    def __init__(self, remover, selector):
        """
        Initialize a :class:`.RemoveMolecules` instance.

        Parameters
        ----------
        remover : :class:`.Selector`
            Selects batches molecules, any molecule selected cannot
            be selected by `selector`.

        selector : :class:`.Selector`
            Selects batches of molecules, not containing any molecules
            selected by `remover`.

        """

        self._remover = remover
        self._selector = selector

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None,
    ):
        remover_batches = self._remover.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
        removed = {
            record
            for batch in remover_batches
            for record in batch
        }
        valid_population = tuple(
            record for record in population if record not in removed
        )
        yield from self._selector.select(
            population=valid_population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
