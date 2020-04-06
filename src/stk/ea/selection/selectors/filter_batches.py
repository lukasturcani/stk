"""
Filter Batches
==============

"""

from .selector import Selector


class FilterBatches(Selector):
    """
    Allows a :class:`.Selector` to select only some batches.

    Examples
    --------
    *Using a Selection Algorithm on a Subset of Batches*

    You only want the :class:`.Best` 10 batches to participate in
    :class:`.Roulette`

    .. code-block:: python

        import stk

        # Select the 10 best batches first, and then run roulette on
        # those.
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

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None,
    ):
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
