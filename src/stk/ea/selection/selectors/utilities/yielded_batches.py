"""
Yielded Batches
===============

"""


class YieldedBatches:
    """
    Keeps track of batches yielded by :meth:`.Selector.select`.

    """

    __slots__ = ('_records', '_batches', '_num', '_key_maker')

    def __init__(self, key_maker):
        self._key_maker = key_maker

        # Has all records yielded by select().
        self._records = set()
        # Has the identity_key() of all batches yielded by select().
        self._batches = set()
        # Counts the total number of times select() has yielded.
        self._num = 0

    def update(self, batch):
        """
        Update tracked data with a new `batch`.

        Parameters
        ----------
        batch : :class:`.Batch`
            A batch yielded by :meth:`.Selector.select`.

        Returns
        -------
        :class:`YieldedBatches`
            The data tracker.

        """

        self._records.update(map(self._key_maker.get_key, batch))
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
