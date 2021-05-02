"""
Random Crosser
==============

"""


import numpy as np


class RandomCrosser:
    """
    Use some other crosser at random.

    Examples
    --------
    *Use One of Several Crossers at Random*

    .. testcode:: use-one-of-several-crossers-at-random

        import stk

        def get_functional_group_type(building_block):
            fg, = building_block.get_functional_groups(0)
            return type(fg)

        def get_num_functional_groups(building_block):
           return building_block.get_num_functional_groups()

        crosser = stk.RandomCrosser(
            crossers=(
                stk.GeneticRecombination(get_functional_group_type),
                stk.GeneticRecombination(get_num_functional_groups),
            ),
            random_seed=12,
        )

        record1 = stk.MoleculeRecord(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='O=CC(C=O)CC=O',
                        functional_groups=[stk.AldehydeFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='NCCN',
                        functional_groups=[stk.PrimaryAminoFactory()],
                    ),
                ),
            ),
        )
        record2 = stk.MoleculeRecord(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='O=CNC(C=O)CC=O',
                        functional_groups=[stk.AldehydeFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='NCNCN',
                        functional_groups=[stk.PrimaryAminoFactory()],
                    ),
                ),
            ),
        )

        # Use one of the component crossers at random.
        crossover_records1 = tuple(crosser.cross((record1, record2)))
        # A different crosser may get selected at random the second,
        # third, etc, time.
        crossover_records2 = tuple(crosser.cross((record1, record2)))

    .. testcode:: use-one-of-several-crossers-at-random
        :hide:

        _expected_results = (
            stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='O=CNC(C=O)CC=O',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='NCNCN',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                    ),
                ),
            ),
            stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='O=CC(C=O)CC=O',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='NCCN',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                    ),
                ),
            ),
            stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='O=CNC(C=O)CC=O',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='NCCN',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                    ),
                ),
            ),
            stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='O=CC(C=O)CC=O',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='NCNCN',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                    ),
                ),
            ),
        )

        def _get_smiles(item):
            if isinstance(item, stk.ConstructedMolecule):
               return stk.Smiles().get_key(item)
            return stk.Smiles().get_key(
                molecule=item.get_molecule_record().get_molecule(),
             )

        _expected_smiles = set(map(_get_smiles, _expected_results))

        _crossover_smiles1 = set(map(_get_smiles, crossover_records1))
        _crossover_smiles2 = set(map(_get_smiles, crossover_records2))
        assert _expected_smiles == _crossover_smiles1
        assert _expected_smiles == _crossover_smiles2

    """

    def __init__(self, crossers, weights=None, random_seed=None):
        """
        Initialize a :class:`.RandomCrosser` instance.

        Parameters
        ----------
        crossers : :class:`tuple`
            Holds instances which have :meth:`.cross` method. The
            :meth:`.cross` method must return an instance of
            :class:`.CrossoverRecord`.

        weights : :class:`tuple` of :class:`float`, optional
            For each crosser, the probability that it will be chosen
            whenever :meth:`.cross` is called.
            If ``None`` all `crossers` will have equal chance of being
            selected.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._crossers = crossers
        self._weights = weights
        self._generator = np.random.RandomState(random_seed)

    def cross(self, records):
        """
        Cross `records`.

        Parameters
        ----------
        records : :class:`iterable` of :class:`.MoleculeRecord`
            The molecule records on which a crossover operation is
            performed.

        Yields
        -------
        :class:`.CrossoverRecord`
            A record of a crossover operation.

        """

        crosser = self._generator.choice(
            a=self._crossers,
            p=self._weights,
        )
        return crosser.cross(records)
