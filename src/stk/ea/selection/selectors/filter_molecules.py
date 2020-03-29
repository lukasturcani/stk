from .selector import Selector


class FilterMolecules(Selector):
    """
    Allows a :class:`.Selector` to select only some molecules.

    Examples
    --------
    You want to use :class:`.Roulette` on the molecules which belong
    to the :class:`.Best` 5 batches of size 3

    .. code-block:: python

        import stk

        population = stk.Population(...)
        selector = stk.FilterMolecules(
            filter=stk.Best(num_batches=5, batch_size=3),
            selector=stk.Roulette(num_batches=20, batch_size=3),
        )
        for batch in selector.select(population):
            # Do stuff with batch. All the molecules in the batch
            # belong to the top 5 batches of size 3. The batch
            # was selected using roulette selection.

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

    def _select(self, population, included_batches, excluded_batches):
        valid_batches = self._filter.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
        valid_pop = population.__class__(*(
            mol for batch in valid_batches for mol in batch
        ))
        valid_pop.set_fitness_values_from_dict(
            fitness_values=population.get_fitness_values(),
        )
        yield from self._selector.select(
            population=valid_pop,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )


