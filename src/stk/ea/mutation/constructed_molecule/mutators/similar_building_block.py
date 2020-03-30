"""
Similar Building Block
======================

"""

import numpy as np

from .mutator import ConstructedMoleculeMutator
from ..record import ConstructedMoleculeMutationRecord
from ....molecule_records import ConstructedMoleculeRecord


class SimilarBuildingBlock(ConstructedMoleculeMutator):
    """
    Substitutes similar building blocks.

    This mutator takes a :class:`.ConstructedMolecule` and substitutes
    the building blocks with the most similar one from a given set.

    Examples
    --------
    *Constructed Molecule Mutation*

    .. code-block:: python

        import stk

        # Create a molecule which is to be mutated.
        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        polymer = stk.ConstructedMoleculeRecord(
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
            return fg.__class__ is stk.PrimaryAmino

        similar_bb = stk.SimilarBuildingBlock(
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

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """

        self._building_blocks = building_blocks
        self._is_replaceable = is_replaceable
        self._name = name
        self._generator = np.random.RandomState(random_seed)

    def _mutate(self, mol):
        """
        Return a mutant of `mol`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule to be mutated.

        Returns
        -------
        mol : :class:`.ConstructedMolecule`
            The mutant.

        """

        if mol not in self._similar_bbs:
            # Maps the mol to a dict. The dict maps each
            # building block of the mol to an iterator.
            # The iterators yield the next most similar molecules in
            # `building_blocks` to the building block.
            self._similar_bbs[mol] = {}

        # Choose the building block which undergoes mutation.
        valid_bbs = [
            bb for bb in mol.building_block_vertices if self._key(bb)
        ]
        chosen_bb = self._generator.choice(valid_bbs)

        # If the building block has not been chosen before, create an
        # iterator yielding similar molecules from `building_blocks`
        # for it.
        if chosen_bb not in self._similar_bbs[mol]:
            self._similar_bbs[mol][chosen_bb] = iter(sorted(
                self._building_blocks,
                key=lambda m: dice_similarity(m, chosen_bb),
                reverse=True
            ))

        try:
            new_bb = next(self._similar_bbs[mol][chosen_bb])
        except StopIteration:
            self._similar_bbs[mol][chosen_bb] = iter(sorted(
                self._building_blocks,
                key=lambda m: dice_similarity(m, chosen_bb),
                reverse=True
            ))
            new_bb = next(self._similar_bbs[mol][chosen_bb])

        # If the most similar molecule in `mols` is itself, then take
        # the next most similar one.
        if new_bb is chosen_bb:
            try:
                new_bb = next(self._similar_bbs[mol][chosen_bb])
            except StopIteration:
                self._similar_bbs[mol][chosen_bb] = iter(sorted(
                    self._building_blocks,
                    key=lambda m: dice_similarity(m, chosen_bb),
                    reverse=True
                ))
                new_bb = next(self._similar_bbs[mol][chosen_bb])

        # Build the new ConstructedMolecule.
        new_bbs = [
            bb for bb in mol.building_block_vertices
            if bb is not chosen_bb
        ]
        new_bbs.append(new_bb)
        return mol.__class__(
            building_blocks=new_bbs,
            topology_graph=mol.topology_graph,
            use_cache=self._use_cache
        )


