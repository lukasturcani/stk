"""
Random Building Block
=====================

"""

import numpy as np

from ....molecule_records import MoleculeRecord
from ...records import MutationRecord
from .mutator import MoleculeMutator


class RandomBuildingBlock(MoleculeMutator):
    """
    Substitutes random building blocks.

    This mutator takes a :class:`.ConstructedMolecule` and substitutes
    the building blocks with one chosen at random from a given set.

    Examples
    --------
    *Constructed Molecule Mutation*

    .. testcode:: constructed-molecule-mutation

        import stk

        # Create a molecule which is to be mutated.
        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        polymer = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear((bb1, bb2), 'AB', 3),
        )

        # Create molecules used to substitute building blocks.
        building_blocks = (
            stk.BuildingBlock(
                smiles='NC[Si]CCN',
                functional_groups=[stk.PrimaryAminoFactory()],
            ),
            stk.BuildingBlock(
                smiles='NCCCCCCCN',
                functional_groups=[stk.PrimaryAminoFactory()],
            ),
            stk.BuildingBlock(
                smiles='NC1CCCCC1N',
                functional_groups=[stk.PrimaryAminoFactory()],
            ),
        )

        # Create the mutator.

        def has_primary_amino_group(building_block):
            fg, = building_block.get_functional_groups(0)
            return type(fg) is stk.PrimaryAmino

        random_bb = stk.RandomBuildingBlock(
            building_blocks=building_blocks,
            is_replaceable=has_primary_amino_group,
        )

        # Mutate a molecule.
        mutation_record1 = random_bb.mutate(polymer)

        # Mutate the molecule a second time.
        mutation_record2 = random_bb.mutate(polymer)

    """

    def __init__(
        self,
        building_blocks,
        is_replaceable,
        name='RandomBuildingBlock',
        random_seed=None,
    ):
        """
        Initialize a :class:`.RandomBuildingBlock` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            A group of molecules which are used to replace building
            blocks in molecules being mutated.

        is_replaceable : :class:`callable`
            A function which takes a :class:`.BuildingBlock` and
            returns ``True`` or ``False``. This function is applied to
            every building block in the molecule being mutated.
            Building blocks which returned ``True`` are liable for
            substitution by one of the molecules in `building_blocks`.

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._building_blocks = building_blocks
        self._is_replaceable = is_replaceable
        self._name = name
        self._generator = np.random.RandomState(random_seed)

    def mutate(self, record):
        # Choose the building block which undergoes mutation.
        replaceable_building_blocks = tuple(filter(
            self._is_replaceable,
            record.get_molecule().get_building_blocks(),
        ))
        replaced_building_block = self._generator.choice(
            a=replaceable_building_blocks,
        )
        # Choose a replacement building block.
        replacement = self._generator.choice(self._building_blocks)

        # Build the new ConstructedMolecule.
        graph = record.get_topology_graph().with_building_blocks({
            replaced_building_block: replacement,
        })
        return MutationRecord(
            molecule_record=MoleculeRecord(graph),
            mutator_name=self._name,
        )
