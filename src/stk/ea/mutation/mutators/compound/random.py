"""
Random Mutator
==============

"""


import numpy as np


class RandomMutator:
    """
    Use some other mutator at random.

    Examples
    --------
    *Use One of Several Mutators at Random*

    .. testcode:: use-one-of-several-mutators-at-random

        import stk

        random_seed = 12
        mutator = stk.RandomMutator(
            mutators=(
                stk.RandomBuildingBlock(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    is_replaceable=lambda building_block: True,
                    random_seed=random_seed,
                ),
                stk.SimilarBuildingBlock(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCCCCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    is_replaceable=lambda building_block: True,
                    random_seed=random_seed,
                ),
                stk.RandomBuildingBlock(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCNCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCNCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    is_replaceable=lambda building_block: True,
                    random_seed=random_seed,
                ),
            ),
        )
        building_block = stk.BuildingBlock(
            smiles='BrCNNCBr',
            functional_groups=[stk.BromoFactory()],
        )
        record = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear(
                building_blocks=(building_block, ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        # Use one of the component mutators at random.
        mutation_record1 = mutator.mutate(record)
        # A different mutator may get selected at random the second,
        # third, etc, time.
        mutation_record2 = mutator.mutate(record)

    .. testcode:: use-one-of-several-mutators-at-random
        :hide:

        _smiles = stk.Smiles()
        assert _smiles.get_key(record.get_molecule()) == 'BrCNNCCNNCBr'

        _molecule1 = (
            mutation_record1.get_molecule_record().get_molecule()
        )
        assert _smiles.get_key(_molecule1) != 'BrCNNCNNCBr'

        _molecule2 = (
            mutation_record2.get_molecule_record().get_molecule()
        )
        assert _smiles.get_key(_molecule2) != 'BrCNNCNNCBr'


    """

    def __init__(self, mutators, weights=None, random_seed=None):
        """
        Initialize a :class:`.RandomMutator` instance.

        Parameters
        ----------
        mutators : :class:`tuple`
            Holds instances which have :meth:`.mutate` method. The
            :meth:`.mutate` method must return an instance of
            :class:`.MutationRecord`.

        weights : :class:`tuple` of :class:`float`, optional
            For each mutator, the probability that it will be chosen
            whenever :meth:`.mutate` is called.
            If ``None`` all `mutators` will have equal chance of being
            selected.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._mutators = mutators
        self._weights = weights
        self._generator = np.random.RandomState(random_seed)

    def mutate(self, record):
        """
        Return a mutant of `record`.

        Parameters
        ----------
        record : :class:`.MoleculeRecord`
            The molecule to be mutated.

        Returns
        -------
        :class:`.MutationRecord`
            A record of the mutation. The exact subclass of
            :class:`.MutationRecord` depends on which mutator was
            used.

        None : :class:`NoneType`
            If `record` cannot be mutated.

        """

        mutator = self._generator.choice(
            a=self._mutators,
            p=self._weights,
        )
        return mutator.mutate(record)
