"""
Selector
========

"""


class Selector:
    """
    An abstract base class for selectors.

    Selectors select batches of molecules from a population.
    Each batch is selected based on its fitness. The fitness of a
    batch is the sum of all fitness values of the molecules in the
    batch. Batches may be of size 1.

    """

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None
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
        :class:`Batch` of :class:`.Molecule`
            A batch of selected molecules.

        """

        # This method can be used to decorate _select in the future.
        yield from self._select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )

    def _select(self, population, included_batches, excluded_batches):
        """
        Select batches of molecules from `population`.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            A collection of molecules from which batches are selected.

        included_batches : :class:`set`
            The identity keys of batches which are allowed to be
            yielded, if ``None`` all batches can be yielded. If not
            ``None`` only batches `included_batches` will be yielded.

        excluded_batches : class:`set`
            The identity keys of batches which are not allowed to be
            yielded. If ``None``, no batch is forbidden from being
            yielded.

        Yields
        ------
        :class:`Batch` of :class:`.Molecule`
            A batch of selected molecules.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()
