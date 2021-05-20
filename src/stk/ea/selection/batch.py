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

    .. testcode:: sorting-batches-by-fitness-value

        import stk

        record1 = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        record2 = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        record3 = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )

        batches = (
            stk.Batch(
                records=(record1, ),
                fitness_values={record1: 1},
                key_maker=stk.Inchi(),
            ),
            stk.Batch(
                records=(record2, ),
                fitness_values={record2: 2},
                key_maker=stk.Inchi(),
            ),
            stk.Batch(
                records=(record3, ),
                fitness_values={record3: 3},
                key_maker=stk.Inchi(),
            ),
        )
        sorted_batches = sorted(batches)

    .. testcode:: sorting-batches-by-fitness-value
        :hide:

        assert (
            sorted_batches[0] < sorted_batches[1] < sorted_batches[2]
        )

    *Comparing Batches by Fitness Value*

    Comparison is also based on fitness value

    .. testcode:: comparing-batches-by-fitness-value

        import stk

        record1 = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        batch1 = stk.Batch(
            records=(record1, ),
            fitness_values={record1: 1},
            key_maker=stk.Inchi(),
        )

        record2 = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        batch2 = stk.Batch(
            records=(record2, ),
            fitness_values={record2: 2},
            key_maker=stk.Inchi(),
        )

        if batch1 < batch2:
            print('batch1 has a smaller fitness value than batch2.')

    .. testoutput:: comparing-batches-by-fitness-value
        :hide:

        batch1 has a smaller fitness value than batch2.

    *Iterating Through Molecule Records in a Batch*

    Batches can be iterated through to get the molecule records in the
    batch

    .. testcode:: iterating-through-molecule-records-in-a-batch

        import stk

        record1 = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        record2 = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        batch = stk.Batch(
            records=(record1, record2),
            fitness_values={record1: 1, record2: 2},
            key_maker=stk.Inchi(),
        )
        for record in batch:
            # Do stuff with record.
            pass

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
