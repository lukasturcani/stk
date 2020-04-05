"""
Selector
========

#. :class:`.AboveAverage`
#. :class:`.Best`
#. :class:`.FilterBatches`
#. :class:`.FilterMoleculeRecords`
#. :class:`.RemoveBatches`
#. :class:`.RemoveMoleculeRecords`
#. :class:`.Roulette`
#. :class:`.StochasticUniversalSampling`
#. :class:`.Tournament`
#. :class:`.Worst`

"""


class Selector:
    """
    An abstract base class for selectors.

    Selectors select batches of molecules from a population.
    Each batch is selected based on its fitness. The fitness of a
    batch is the sum of all fitness values of the molecules in the
    batch. Batches may be of size 1.

    See Also
    --------
    :class:`.Batch`

    Examples
    --------
    *Subclass Implementation*

    The source codeof the classes listed in :mod:`.selector` can serve
    as good examples.

    """

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None,
    ):
        """
        Select batches of molecules from `population`.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            A collection of molecules from which batches are selected.

        included_batches : :class:`set`, optional
            The identity keys of batches which are allowed to be
            yielded, if ``None`` all batches can be yielded. If not
            ``None`` only batches `included_batches` will be yielded.

        excluded_batches : class:`set`, optional
            The identity keys of batches which are not allowed to be
            yielded. If ``None``, no batch is forbidden from being
            yielded.

        Yields
        ------
        :class:`Batch` of :class:`.MoleculeRecord`
            A batch of selected molecule records.

        """

        raise NotImplementedError()
