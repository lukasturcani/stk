
class FilterBatches(Selector):
    """
    Allows a :class:`.Selector` to select only some batches.

    Examples
    --------
    You only want the :class:`.Best` 10 batches to participate in
    :class:`.Roulette`

    .. code-block:: python

        import stk

        population = stk.Population(...)
        selector = stk.FilterBatches(
            filter=stk.Best(10),
            selector=stk.Roulette(7),
        )
        for batch in selector.select(population):
            # Do stuff with batch. It is one of the 10 best batches and
            # was selected using roulette selection.

    """

    def __init__(self, filter, selector):
        """
        Initialize a :class:`.FilterBatches` instance.

        Parameters
        ----------
        filter : :class:`.Selector`
            Selects batches which can be yielded by `selector`.

        selector : :class:`.Selector`
            Selects batches, but only if they were also selected by
            `filter`.

        """

        self._filter = filter
        self._selector = selector

    def _select(self, population, included_batches, excluded_batches):
        allowed = {
            batch.get_identity_key()
            for batch in self._filter.select(
                population=population,
                included_batches=included_batches,
                excluded_batches=excluded_batches,
            )
        }
        yield from self._selector.select(
            population=population,
            included_batches=allowed,
            excluded_batches=excluded_batches,
        )


