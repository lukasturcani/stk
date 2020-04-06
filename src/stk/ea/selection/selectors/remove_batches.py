"""
Remove Batches
==============

"""

from .selector import Selector


class RemoveBatches(Selector):
    """
    Prevents a :class:`.Selector` from selecting some batches.

    Examples
    --------
    *Removing Batches From Selection*

    You want to use :class:`.Roulette` selection on all but the
    5 :class:`.Worst` batches

    .. code-block:: python

        import stk

        selector = stk.RemoveBatches(
            remover=stk.Worst(5),
            selector=stk.Roulette(20),
        )
        for batch in selector.select(population):
            # Do stuff with batch. It was selected with roulette
            # selection and is not one of the worst 5 batches.

    """

    def __init__(self, remover, selector):
        """
        Initialize a :class:`.RemoveBatches` instance.

        Parameters
        ----------
        remover : :class:`.Selector`
            Selects batches of molecules, which cannot be yielded by
            `selector`.

        selector : :class:`.Selector`
            Selects batches of molecules, except those selected by
            `remover`.

        """

        self._remover = remover
        self._selector = selector

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None,
    ):
        removed_batches = {
            batch.get_identity_key()
            for batch in self._remover.select(
                population=population,
                included_batches=included_batches,
                excluded_batches=excluded_batches,
            )
        }
        if excluded_batches is not None:
            removed_batches |= excluded_batches

        yield from self._selector.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=removed_batches,
        )
