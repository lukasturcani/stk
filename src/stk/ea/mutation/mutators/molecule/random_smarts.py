"""
Random Functional Group
=====================

"""

import numpy as np


from .mutator import MoleculeMutator
from ...records import MutationRecord
from ....molecule_records import MoleculeRecord
from rdkit import Chem
from rdkit.Chem.rdchem import BondType
from .....molecular import BuildingBlock


class RandomSmartsFunctionalGroup(MoleculeMutator):
    """
    Substitutes functional groups within building blocks.

    This mutator takes a :class:`ConstructedMolecule` and substitutes
    the existing building blocks with ones containing
    new functional groups from a given set of
    `class`:`SmartsFunctionalGroupFactory`.
    This function uses the bonder ids in the provided functional
    groups to attach to the original molecule.
    """
    def __init__(
        self,
        functional_groups,
        is_replaceable_fg,
        is_replaceable_bb,
        replacement_count='all',
        name='RandomFunctionalGroup',
        random_seed=None,
    ):
        """
        Initialize a :class:`.RandomFunctionalGroup` instance.

        Parameters
        ----------
        functional_groups : :class:`tuple` of
        :class`FunctionalGroupFactories`
            A group of :class:`SmartsFunctionalGroupsFactory`
            which are used to build :class:`FunctionalGroup` that
            replace :class:`FunctionalGroup` in the existing
            building blocks.

        is_replaceable_bb : :class:`callable`
            A function which takes a :class:`.BuildingBlock` and
            returns ``True`` or ``False``. This function is applied to
            every building block in the molecule being mutated.
            Building blocks which returned ``True`` are liable for
            substitution by one of the molecules in `building_blocks`.

        is_replaceable_fg : :class:`callable`
            A function which takes a :class:``
             and returns ``True`` or ``False``.
            This function is applied to
            every building block in the molecule being mutated.
            Functional groups in the building block which returned
            ``True`` are liable for
            substitution by one of the functional groups provided in
             the `functional groups` argument.

        replacement_count : :class:`str` or :class:`int`, optional
            Specifies the number of functional groups in a building
            block to relpace.
            If `all`, all functional groups within the building block
            to replace.

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """
        self._functional_groups = functional_groups
        self._is_replaceable_bb = is_replaceable_bb
        self._is_replaceable_fg = is_replaceable_fg
        self._name = name
        self._replacement_count = replacement_count
        self._generator = np.random.RandomState(random_seed)

    def mutate(self, record):
        replaceable_building_blocks = tuple(filter(
            self._is_replaceable_bb,
            record.get_molecule().get_building_blocks(),
        ))
        replaced_building_block = self._generator.choice(
            a=replaceable_building_blocks,
        )
        replaced_fgs = tuple(filter(
            self._is_replaceable_fg,
            replaced_building_block.get_functional_groups()
            )
        )
        # Creates tuple of replaceable functional groups
        if isinstance(self._replacement_count, int):
            replaced_fgs = self._generator.choice(
                a=replaced_fgs,
                size=self._replacement_count
            )
        elif self._replacement_count != 'all':
            raise TypeError('Unsupported type for replacement count.')
        # Select a new factory to replace all
        # selected functional groups
        fg_replacer = self._generator.choice(
            a=self._functional_groups
        )

        rd_mol = replaced_building_block.to_rdkit_mol()
        # Flat list of bonder ids to keep
        # Bonder ids are from functional groups to be replaced
        replaced_bonder_ids = [
            bonder_id
            for fg in replaced_fgs
            for bonder_id in fg.get_bonder_ids()
        ]
        # Flat list of atom ids to remove
        # Atom ids are from the original functional groups
        deleted_ids = [
            atom_id
            for fg in replaced_fgs
            for atom_id in fg.get_atom_ids()
            if atom_id not in replaced_bonder_ids
        ]
        # Perform molecule editting
        editable_mol = Chem.rdchem.RWMol(rd_mol)
        for atom_id in sorted(deleted_ids, reverse=True):
            editable_mol.RemoveAtom(atom_id)
        # return editable_mol
        new_fg = set_mapped_atoms(
            Chem.MolFromSmarts(fg_replacer._smarts)
        )
        if isinstance(fg_replacer._bonders, int):
            fg_bonder_ids = [fg_replacer._bonders]
        else:
            fg_bonder_ids = [
                bonder_id for bonder_id in fg_replacer._bonders
            ]
        # Insert new functional groups
        for fg in [new_fg]*len(replaced_fgs):
            editable_mol.InsertMol(fg)
        # Get inserted functional group's bonder id
        new_bonder_ids = [
            atom.GetIdx()
            for atom in
            get_mapped_atoms(editable_mol, fg_bonder_ids)
        ]
        print(fg_bonder_ids)
        print(new_bonder_ids)
        print(replaced_bonder_ids)
        for a1_id, a2_id in zip(replaced_bonder_ids, new_bonder_ids):
            print(a1_id, a2_id)
            editable_mol.AddBond(a1_id, a2_id, BondType.SINGLE)
        # Build new building block
        replacement = BuildingBlock.init_from_rdkit_mol(
            editable_mol,
            functional_groups=[fg_replacer]
        )
        return (new_fg, replacement, replaced_fgs)
        graph = record.get_topology_graph().with_building_blocks({
            replaced_building_block: replacement
        })
        return MutationRecord(
            molecule_record=MoleculeRecord(graph),
            mutator_name=self._name,
        )


def get_mapped_atoms(mol, mapper_ids):
    """
    Filters for atoms in a :class: `RDKit.Mol` with a specific mapping
    number.
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() in mapper_ids:
            yield atom


def set_mapped_atoms(mol):
    """
    Sets the atoms map number for all atoms in a molecule with
    the default atom id for that atom.
    """
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol
