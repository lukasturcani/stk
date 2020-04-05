"""

"""


class YieldedBatches:
    """
    Keeps track of batches yielded by :meth:`.Selector._select`.

    Each :meth:`.Selector._select` call should be paired with a
    new :class:`._YieldedData` instance, which should be updated
    each time a new batch is yielded.

    """

    __slots__ = ['_mols', '_batches', '_num']

    def __init__(self):
        # Has all molecules yielded by _select().
        self._mols = set()
        # Has the identity_key() of all batches yielded by select().
        self._batches = set()
        # Counts the total number of times _select() has yielded.
        self._num = 0

    def update(self, batch):
        """
        Update tracked data with a new `batch`.

        Parameters
        ----------
        batch : :class:`.Batch`
            A batch yielded by :meth:`.Selector._select`.

        Returns
        -------
        :class:`_YieldedData`
            The data tracker.

        """

        self._mols.update(batch)
        self._batches.add(batch.get_identity_key())
        self._num += 1
        return self

    def get_num(self):
        """
        Get the number of times :meth:`.Selector._select` has yielded.

        Returns
        -------
        :class:`int`
            The total number of times :meth:`.Selector._select` has
            yielded.

        """

        return self._num

    def is_yielded_batch(self, batch):
        """
        Check if `batch` has already been yielded.

        Parameters
        ----------
        batch : :class:`.Batch`
            The batch to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `batch` has already been yielded.

        """

        return batch.get_identity_key() in self._batches

    def is_unyielded_batch(self, batch):
        """
        Check if `batch` has not been yielded.

        Parameters
        ----------
        batch : :class:`.Batch`
            The batch to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `batch` has not been yielded.

        """

        return batch.get_identity_key() not in self._batches

    def has_yielded_mols(self, batch):
        """
        Check if `batch` contains any previously yielded molecules.

        Parameters
        ----------
        batch : :class:`.Batch`
            The batch to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `batch` contains any molecules which have
            previously been yielded.

        """

        return any(mol in self._mols for mol in batch)

    def has_no_yielded_mols(self, batch):
        """
        Check if `batch` consists only of unyielded molecules.

        Parameters
        ----------
        batch : :class:`.Batch`
            The batch to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `batch` does not have any previously yielded
            molecules.

        """

        return all(mol not in self._mols for mol in batch)
