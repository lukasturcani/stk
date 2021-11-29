"""
Similar Building Block
======================

"""

from functools import partial

import numpy as np

from stk.molecular import Inchi
from stk.utilities import dice_similarity

from ....molecule_records import MoleculeRecord
from ...records import MutationRecord
from .mutator import MoleculeMutator


class SimilarBuildingBlock(MoleculeMutator):
    """
    Substitutes similar building blocks.

    This mutator takes a :class:`.ConstructedMolecule` and substitutes
    the building blocks with the most similar one from a given set.
    Repeated mutations on the same molecule will substituted the next
    most similar molecule from the set.

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

        similar_bb = stk.SimilarBuildingBlock(
            building_blocks=building_blocks,
            is_replaceable=has_primary_amino_group,
        )

        # Mutate a molecule.
        mutation_record1 = similar_bb.mutate(polymer)

        # Mutate the molecule a second time.
        mutation_record2 = similar_bb.mutate(polymer)

    """

    def __init__(
        self,
        building_blocks,
        is_replaceable,
        key_maker=Inchi(),
        name='SimilarBuildingBlock',
        random_seed=None,
    ):
        """
        Initialize a :class:`.SimilarBuildingBlock` instance.

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

        key_maker : :class:`.MoleculeKeyMaker`, optional
            Molecules which return the same key, will iterate through
            the same set of similar molecules.

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """

        self._building_blocks = building_blocks
        self._is_replaceable = is_replaceable
        self._key_maker = key_maker
        self._name = name
        self._generator = np.random.RandomState(random_seed)
        self._similar_building_blocks = {}

    def mutate(self, record):
        key = self._key_maker.get_key(record.get_molecule())
        if key not in self._similar_building_blocks:
            # Maps the key to a dict. The dict maps each
            # building block to an iterator.
            # The iterators yield the next most similar molecules in
            # `building_blocks` to the building block.
            self._similar_building_blocks[key] = {}

        similar_building_blocks = self._similar_building_blocks[key]

        # Choose the building block which undergoes mutation.
        replaceable_building_blocks = tuple(filter(
            self._is_replaceable,
            record.get_molecule().get_building_blocks(),
        ))
        replaced_building_block = self._generator.choice(
            a=replaceable_building_blocks,
        )

        # If the building block has not been chosen before, create an
        # iterator yielding similar molecules from `building_blocks`
        # for it.
        replaced_key = self._key_maker.get_key(replaced_building_block)
        if replaced_key not in similar_building_blocks:
            similar_building_blocks[replaced_key] = iter(sorted(
                self._building_blocks,
                key=partial(dice_similarity, replaced_building_block),
                reverse=True,
            ))

        try:
            replacement = next(similar_building_blocks[replaced_key])
        except StopIteration:
            similar_building_blocks[replaced_key] = iter(sorted(
                self._building_blocks,
                key=partial(dice_similarity, replaced_building_block),
                reverse=True,
            ))
            replacement = next(similar_building_blocks[replaced_key])

        # If the most similar molecule in `building_blocks` is itself,
        # then take the next most similar one.
        if self._key_maker.get_key(replacement) == replaced_key:
            try:
                replacement = next(
                    similar_building_blocks[replaced_key]
                )
            except StopIteration:
                similar_building_blocks[replaced_key] = iter(sorted(
                    self._building_blocks,
                    key=partial(
                        dice_similarity,
                        replaced_building_block,
                    ),
                    reverse=True,
                ))
                replacement = next(
                    similar_building_blocks[replaced_key]
                )

        # Build the new ConstructedMolecule.
        graph = record.get_topology_graph().with_building_blocks({
            replaced_building_block: replacement,
        })
        return MutationRecord(
            molecule_record=MoleculeRecord(graph),
            mutator_name=self._name,
        )
