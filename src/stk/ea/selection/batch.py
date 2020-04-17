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

    __slots__ = ('_records', '_fitness_value', '_identity_key')

    def __init__(self, records, fitness_values, key_maker):
        """
        Initialize a :class:`.Batch`.

        Parameters
        ----------
        records : :class:`tuple` of :class:`.MoleculeRecord`
            The molecule records which are part of the batch.

        fitness_values : :class:`dict`
            Maps each :class:`.MoleculeRecord` in `records` to the
            fitness value which should be used for it.

        key_maker : :class:`.MoleculeKeyMaker`
            Used to make keys for molecules, which are used to
            determine the identity key of the batch. If two
            batches have the same molecule keys, the same number of
            times, they will have the same identity key.

        """

        self._records = records
        self._fitness_value = sum(map(fitness_values.get, records))
        molecules = (record.get_molecule() for record in records)
        self._identity_key = frozenset(
            Counter(map(key_maker.get_key, molecules)).items()
        )

    def get_size(self):
        """
        Get the number of molecules in the batch.

        Returns
        -------
        :class:`int`
            The number of molecules in the batch.

        """

        return len(self._records)

    def get_fitness_value(self):
        """
        Get the fitness value of the batch.

        Returns
        -------
        :class:`float`
            The fitness value.

        """

        return self._fitness_value

    def get_identity_key(self):
        """
        Get the identity key of the batch.

        If two batches hold the same molecules, the same number of
        times, they will have the same identity key.

        Returns
        -------
        :class:`object`
            A hashable object which can be used to compare if two
            batches have the same identity.

        """

        return self._identity_key

    def __iter__(self):
        return iter(self._records)

    def __getitem__(self, index):
        return self._records[index]

    def __eq__(self, other):
        return self._fitness_value == other._fitness_value

    def __gt__(self, other):
        return self._fitness_value > other._fitness_value

    def __ge__(self, other):
        return self._fitness_value >= other._fitness_value

    def __lt__(self, other):
        return self._fitness_value < other._fitness_value

    def __le__(self, other):
        return self._fitness_value <= other._fitness_value

    def __repr__(self):
        return f'Batch({self._fitness_value})'

    def __str__(self):
        return repr(self)
