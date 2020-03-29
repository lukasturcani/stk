from .selector import Selector


class RemoveMolecules(Selector):
    """
    Prevents a :class:`.Selector` from selecting some molecules.

    Examples
    --------
    You want to prevent any of the molecules in the :class:`.Best`
    5 batches from being selected by :class:`.Roulette`.

    .. code-block:: python

        import stk

        population = stk.Population(...)
        selector = stk.RemoveMolecules(
            remover=stk.Best(num_batches=5, batch_size=3),
            selector=stk.Roulette(num_batches=20, batch_size=3),
        )

        for batch in selector.select(population):
            # Do stuff with batch. The batch is guaranteed not to
            # contain any molecules which are found in the best 5
            # batches of size 3.

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

    def _select(self, population, included_batches, excluded_batches):
        remover_batches = self._remover.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
        removed = {mol for batch in remover_batches for mol in batch}
        valid_pop = population.__class__(*(
            mol for mol in population if mol not in removed
        ))
        valid_pop.set_fitness_values_from_dict(
            fitness_values=population.get_fitness_values(),
        )
        yield from self._selector.select(
            population=valid_pop,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )


