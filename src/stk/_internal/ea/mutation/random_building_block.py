import typing
from collections.abc import Iterable

import numpy as np

from stk._internal.building_block import BuildingBlock
from stk._internal.ea.molecule_records.molecule import MoleculeRecord
from stk._internal.ea.mutation.record import MutationRecord


class RandomBuildingBlock:
    """
    Substitutes random building blocks.

    This mutator takes a :class:`.ConstructedMolecule` and substitutes
    the building blocks with one chosen at random from a given set.

    Examples:

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
        building_blocks: Iterable[BuildingBlock],
        is_replaceable: typing.Callable[[BuildingBlock], bool],
        name: str = "RandomBuildingBlock",
        random_seed: int | np.random.Generator | None = None,
    ) -> None:
        """
        Parameters:
            building_blocks (list[BuildingBlock]):
                A group of molecules which are used to replace building
                blocks in molecules being mutated.

            is_replaceable:
                This function is applied to every building block in
                the molecule being mutated. Building blocks
                which returned ``True`` are liable for substitution
                by one of the molecules in `building_blocks`.

            name:
                A name to help identify the mutator instance.

            random_seed:
                The random seed to use.
        """
        if random_seed is None or isinstance(random_seed, int):
            random_seed = np.random.default_rng(random_seed)

        self._building_blocks = tuple(building_blocks)
        self._is_replaceable = is_replaceable
        self._name = name
        self._generator = random_seed

    def mutate(self, record: MoleculeRecord) -> MutationRecord[MoleculeRecord]:
        """
        Return a mutant of `record`.

        Parameters:
            record:
                The molecule to be mutated.

        Returns:
            A record of the mutation.
        """
        # Choose the building block which undergoes mutation.
        replaceable_building_blocks = tuple(
            filter(
                self._is_replaceable,
                record.get_molecule().get_building_blocks(),
            )
        )
        replaced_building_block = self._generator.choice(
            a=replaceable_building_blocks,
        )
        # Choose a replacement building block.
        replacement = self._generator.choice(
            self._building_blocks,  # type: ignore
        )

        # Build the new ConstructedMolecule.
        graph = record.get_topology_graph().with_building_blocks(
            {
                replaced_building_block: replacement,
            }
        )
        return MutationRecord(
            molecule_record=MoleculeRecord(graph),
            mutator_name=self._name,
        )
