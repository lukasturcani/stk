"""
Batch
=====

"""

from collections import Counter


class Batch:
    """
    Represents a batch of molecule records.

    Batches can be compared, the comparison is based on their
    fitness values. Batches can also be iterated through, this
    iterates through all the records in the batch.

    Examples
    --------
    *Sorting Batches by Fitness Value*

    Sorting batches causes them to be sorted by fitness value.

    .. code-block:: python

        batches = (Batch(...), Batch(...), Batch(...))
        sorted_batches = sorted(batches)

    *Comparing Batches by Fitness Value*

    Comparison is also based on fitness value

    .. code-block:: python

        batch1 = Batch(...)
        batch2 = Batch(...)
        if batch1 > batch2:
            print('batch1 has a larger fitness value than batch2.')

    *Iterating Through Molecule Records in a Batch*

    Batches can be iterated through to get the molecule records in the
    batch

    .. code-block:: python

        batch = Batch(...)
        for record in batch:
            # Do stuff with record.

    """

    __slots__ = ['_mols', '_fitness', '_identity_key']

    def __init__(self, mols, fitness_values):
        """
        Initialize a :class:`.Batch`.

        Parameters
        ----------
        mols : :class:`tuple` of :class:`.Molecule`
            The molecules which are part of the batch.

        fitness_values : :class:`dict`
            Maps each molecule in `mols` to its fitness value.

        """

        self._mols = mols
        self._fitness = sum(fitness_values[mol] for mol in mols)
        self._identity_key = frozenset(Counter(mols).items())

    def get_size(self):
        """
        Get the number of molecules in the batch.

        Returns
        -------
        :class:`int`
            The number of molecules in the batch.

        """

        return len(self._mols)

    def get_fitness(self):
        """
        Get the fitness value of the batch.

        Returns
        -------
        :class:`float`
            The fitness value.

        """

        return self._fitness

    def get_identity_key(self):
        """
        Get the identity key of the batch.

        If two batches hold the same molecules, the same number of
        times, they will have the same identity key.

        Returns
        -------
        :class:`object`
            A hashable object which can be used to compare if two
            batches are the same.

        """

        return self._identity_key

    def __iter__(self):
        return iter(self._mols)

    def __eq__(self, other):
        return self._fitness == other._fitness

    def __gt__(self, other):
        return self._fitness > other._fitness

    def __ge__(self, other):
        return self._fitness >= other._fitness

    def __lt__(self, other):
        return self._fitness < other._fitness

    def __le__(self, other):
        return self._fitness <= other._fitness

    def __repr__(self):
        return f'Batch({", ".join(str(m) for m in self._mols)})'

    def __str__(self):
        return repr(self)


