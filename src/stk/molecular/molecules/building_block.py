"""
Building Block
==============

"""

from __future__ import annotations

import logging
import os
from collections.abc import Collection
from functools import partial
from typing import Iterable, Iterator, Optional, Union

import numpy as np
import rdkit.Chem.AllChem as rdkit

from ...utilities import flatten, remake
from ..atoms import Atom
from ..bonds import Bond
from ..functional_groups import FunctionalGroup
from .molecule import Molecule

logger = logging.getLogger(__name__)


class BuildingBlock(Molecule):
    """
    Represents a building block of a :class:`.ConstructedMolecule`.

    A :class:`BuildingBlock` can represent either an entire molecule or
    a molecular fragments used to construct a
    :class:`.ConstructedMolecule`. The building block uses
    :class:`.FunctionalGroup` instances to identify which atoms are
    modified during construction.

    """

    # Maps file extensions to functions which can be used to
    # create an rdkit molecule from that file type.
    _init_funcs = {
        '.mol': partial(
            rdkit.MolFromMolFile,
            sanitize=False,
            removeHs=False
        ),

        '.sdf': partial(
            rdkit.MolFromMolFile,
            sanitize=False,
            removeHs=False
        ),

        '.pdb': partial(
            rdkit.MolFromPDBFile,
            sanitize=False,
            removeHs=False,
            proximityBonding=False,
        ),
    }

    def __init__(
        self,
        smiles,
        functional_groups=(),
        placer_ids=None,
        position_matrix=None,
    ):
        """
        Initialize a :class:`.BuildingBlock`.

        Notes
        -----
        The molecule is given 3D coordinates with
        :func:`rdkit.ETKDGv2`.

        Parameters
        ----------
        smiles : :class:`str`
            A SMILES string of the molecule.

        functional_groups : :class:`iterable`, optional
            An :class:`iterable` of :class:`.FunctionalGroup` or
            :class:`.FunctionalGroupFactory` or both.
            :class:`.FunctionalGroup` instances are added to the
            building block and :class:`.FunctionalGroupFactory`
            instances are used to create :class:`.FunctionalGroup`
            instances the building block should hold.
            :class:`.FunctionalGroup` instances are used to identify
            which atoms are modified during
            :class:`.ConstructedMolecule` construction.

        placer_ids : :class:`tuple` of :class:`int`, optional
            The ids of *placer* atoms. These are the atoms which should
            be used for calculating the position of the building block.
            Depending on the values passed to `placer_ids`,
            and the functional groups in the building block, different
            *placer* ids will be used by the building block.

            #. `placer_ids` is passed to the initializer: the passed
               *placer* ids will be used by the building block.

            #. `placer_ids` is ``None`` and the building block has
               functional groups: The *placer* ids of the functional
               groups will be used as the *placer* ids of the building
               block.

            #. `placer_ids` is ``None`` and `functional_groups` is
               empty. All atoms of the molecule will be used for
               *placer* ids.

        position_matrix : :class:`numpy.ndarray`, optional
            The position matrix the building block should use. If
            ``None``, :func:`rdkit.ETKDGv2` will be used to calculate
            it.

        Raises
        ------
        :class:`RuntimeError`
            If embedding the molecule fails.

        """

        molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
        if position_matrix is None:
            params = rdkit.ETKDGv2()
            random_seed = 4
            params.randomSeed = random_seed
            if rdkit.EmbedMolecule(molecule, params) == -1:
                raise RuntimeError(
                    f'Embedding with seed value of {random_seed} '
                    'failed.'
                )
            rdkit.Kekulize(molecule)
        else:
            # Make sure the position matrix always holds floats.
            position_matrix = np.array(
                position_matrix,
                dtype=np.float64,
            )
            conformer = rdkit.Conformer(molecule.GetNumAtoms())
            for atom_id, position in enumerate(position_matrix):
                conformer.SetAtomPosition(atom_id, position)
            molecule.AddConformer(conformer)

        self._init_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
            placer_ids=placer_ids,
        )

    @classmethod
    def init_from_molecule(
        cls,
        molecule,
        functional_groups=(),
        placer_ids=None,
    ):
        """
        Initialize from a :class:`.Molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to initialize from.

        functional_groups : :class:`iterable`, optional
            An :class:`iterable` of :class:`.FunctionalGroup` or
            :class:`.FunctionalGroupFactory` or both.
            :class:`.FunctionalGroup` instances are added to the
            building block and :class:`.FunctionalGroupFactory`
            instances are used to create :class:`.FunctionalGroup`
            instances the building block should hold.
            :class:`.FunctionalGroup` instances are used to identify
            which atoms are modified during
            :class:`.ConstructedMolecule` construction.

        placer_ids : :class:`tuple` of :class:`int`, optional
            The ids of *placer* atoms. These are the atoms which should
            be used for calculating the position of the building block.
            Depending on the values passed to `placer_ids`,
            and the functional groups in the building block, different
            *placer* ids will be used by the building block.

            #. `placer_ids` is passed to the initializer: the passed
               *placer* ids will be used by the building block.

            #. `placer_ids` is ``None`` and the building block has
               functional groups: The *placer* ids of the functional
               groups will be used as the *placer* ids of the building
               block.

            #. `placer_ids` is ``None`` and `functional_groups` is
               empty. All atoms of the molecule will be used for
               *placer* ids.

        Returns
        -------
        :class:`.BuildingBlock`
             The building block. It will have the same atoms, bonds and
             atomic positions as `molecule`.

        """

        return cls.init(
            atoms=tuple(molecule.get_atoms()),
            bonds=tuple(molecule.get_bonds()),
            position_matrix=molecule.get_position_matrix(),
            functional_groups=functional_groups,
            placer_ids=placer_ids,
        )

    @classmethod
    def init_from_vabene_molecule(
        cls,
        molecule,
        functional_groups=(),
        placer_ids=None,
        position_matrix=None,
    ):
        """
        Initialize from a :mod:`vabene.Molecule`.

        Notes
        -----
        The molecule is given 3D coordinates with
        :func:`rdkit.ETKDGv2()`.

        Parameters
        ----------
        molecule : :class:`vabene.Molecule`
            The molecule from which to initialize.

        functional_groups : :class:`iterable`, optional
            An :class:`iterable` of :class:`.FunctionalGroup` or
            :class:`.FunctionalGroupFactory` or both.
            :class:`.FunctionalGroup` instances are added to the
            building block and :class:`.FunctionalGroupFactory`
            instances are used to create :class:`.FunctionalGroup`
            instances the building block should hold.
            :class:`.FunctionalGroup` instances are used to identify
            which atoms are modified during
            :class:`.ConstructedMolecule` construction.

        placer_ids : :class:`tuple` of :class:`int`, optional
            The ids of *placer* atoms. These are the atoms which should
            be used for calculating the position of the building block.
            Depending on the values passed to `placer_ids`,
            and the functional groups in the building block, different
            *placer* ids will be used by the building block.

            #. `placer_ids` is passed to the initializer: the passed
               *placer* ids will be used by the building block.

            #. `placer_ids` is ``None`` and the building block has
               functional groups: The *placer* ids of the functional
               groups will be used as the *placer* ids of the building
               block.

            #. `placer_ids` is ``None`` and `functional_groups` is
               empty. All atoms of the molecule will be used for
               *placer* ids.

        position_matrix : :class:`numpy.ndarray`, optional
            The position matrix the building block should use. If
            ``None``, :func:`rdkit.ETKDGv2` will be used to calculate
            it.

        Returns
        -------
        :class:`.BuildingBlock`
             The building block.

        Raises
        ------
        :class:`RuntimeError`
            If embedding the molecule fails.

        """

        editable = rdkit.EditableMol(rdkit.Mol())
        for atom in molecule.get_atoms():
            rdkit_atom = rdkit.Atom(atom.get_atomic_number())
            rdkit_atom.SetFormalCharge(atom.get_charge())
            editable.AddAtom(rdkit_atom)

        for bond in molecule.get_bonds():
            editable.AddBond(
                beginAtomIdx=bond.get_atom1_id(),
                endAtomIdx=bond.get_atom2_id(),
                order=rdkit.BondType(bond.get_order()),
            )

        rdkit_molecule = editable.GetMol()
        rdkit.SanitizeMol(rdkit_molecule)
        rdkit_molecule = rdkit.AddHs(rdkit_molecule)

        if position_matrix is None:
            params = rdkit.ETKDGv2()
            random_seed = 4
            params.randomSeed = random_seed
            if rdkit.EmbedMolecule(rdkit_molecule, params) == -1:
                raise RuntimeError(
                    f'Embedding with seed value of {random_seed} '
                    'failed.'
                )
        else:
            # Make sure the position matrix always holds floats.
            position_matrix = np.array(
                position_matrix,
                dtype=np.float64,
            )
            conformer = rdkit.Conformer(rdkit_molecule.GetNumAtoms())
            for atom_id, position in enumerate(position_matrix):
                conformer.SetAtomPosition(atom_id, position)
            rdkit_molecule.AddConformer(conformer)

        rdkit.Kekulize(rdkit_molecule)
        return cls.init_from_rdkit_mol(
            molecule=rdkit_molecule,
            functional_groups=functional_groups,
            placer_ids=placer_ids,
        )

    @classmethod
    def init(
        cls,
        atoms,
        bonds,
        position_matrix,
        functional_groups=(),
        placer_ids=None,
    ):
        """
        Initialize a :class:`.BuildingBlock` from its components.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms of the building block.

        bonds : :class:`tuple` of :class:`.Bond`
            The bonds of the building block.

        position_matrix : :class:`numpy.ndarray`
            An ``(n, 3)`` position matrix of the building block.

        functional_groups : :class:`iterable`, optional
            An :class:`iterable` holding the :class:`.FunctionalGroup`
            instances the building block should have, and / or
            :class:`.FunctionalGroupFactory` instances used for
            creating them.

        placer_ids : :class:`tuple` of :class:`int`, optional
            The ids of *placer* atoms. These are the atoms which should
            be used for calculating the position of the building block.
            Depending on the values passed to `placer_ids`,
            and the functional groups in the building block, different
            *placer* ids will be used by the building block.

            #. `placer_ids` is passed to the initializer: the passed
               *placer* ids will be used by the building block.

            #. `placer_ids` is ``None`` and the building block has
               functional groups: The *placer* ids of the functional
               groups will be used as the *placer* ids of the building
               block.

            #. `placer_ids` is ``None`` and `functional_groups` is
               empty. All atoms of the molecule will be used for
               *placer* ids.

        Returns
        -------
        :class:`.BuildingBlock`
            The building block.

        """

        building_block = cls.__new__(cls)
        Molecule.__init__(
            self=building_block,
            atoms=atoms,
            bonds=bonds,
            position_matrix=position_matrix,
        )
        functional_groups = building_block._extract_functional_groups(
            functional_groups=functional_groups,
        )
        building_block._with_functional_groups(functional_groups)
        building_block._placer_ids = (
            building_block._normalize_placer_ids(
                placer_ids=placer_ids,
                functional_groups=building_block._functional_groups,
            )
        )
        building_block._core_ids = frozenset(
            building_block._get_core_ids(
                functional_groups=building_block._functional_groups,
            )
        )
        return building_block

    @classmethod
    def init_from_file(
        cls,
        path,
        functional_groups=(),
        placer_ids=None,
    ):
        """
        Initialize from a file.

        Parameters
        ----------
        path : :class:`str`
            The path to a molecular structure file. Supported file
            types are:

                #. ``.mol``, ``.sdf`` - MDL V3000 MOL file
                #. ``.pdb`` - PDB file

        functional_groups : :class:`iterable`, optional
            An :class:`iterable` of :class:`.FunctionalGroup` or
            :class:`.FunctionalGroupFactory` or both.
            :class:`.FunctionalGroup` instances are added to the
            building block and :class:`.FunctionalGroupFactory`
            instances are used to create :class:`.FunctionalGroup`
            instances the building block should hold.
            :class:`.FunctionalGroup` instances are used to identify
            which atoms are modified during
            :class:`.ConstructedMolecule` construction.

        placer_ids : :class:`tuple` of :class:`int`, optional
            The ids of *placer* atoms. These are the atoms which should
            be used for calculating the position of the building block.
            Depending on the values passed to `placer_ids`,
            and the functional groups in the building block, different
            *placer* ids will be used by the building block.

            #. `placer_ids` is passed to the initializer: the passed
               *placer* ids will be used by the building block.

            #. `placer_ids` is ``None`` and the building block has
               functional groups: The *placer* ids of the functional
               groups will be used as the *placer* ids of the building
               block.

            #. `placer_ids` is ``None`` and `functional_groups` is
               empty. All atoms of the molecule will be used for
               *placer* ids.

        Returns
        -------
        :class:`.BuildingBlock`
            The building block.

        Raises
        ------
        :class:`ValueError`
            If the file type cannot be used for initialization.

        """

        _, extension = os.path.splitext(path)

        if extension not in cls._init_funcs:
            raise ValueError(
                f'Unable to initialize from "{extension}" files.'
            )
        # This remake needs to be here because molecules loaded
        # with rdkit often have issues, because rdkit tries to do
        # bits of structural analysis like stereocenters. Remake
        # gets rid of all this problematic metadata.
        molecule = remake(cls._init_funcs[extension](path))

        return cls.init_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
            placer_ids=placer_ids,
        )

    @classmethod
    def init_from_rdkit_mol(
        cls,
        molecule,
        functional_groups=(),
        placer_ids=None,
    ):
        """
        Initialize from an :mod:`rdkit` molecule.

        Notes
        -----
        For :mod:`rdkit` molecules with non-integer bond orders,
        such as 1.5, the molecule should be kekulized prior to
        calling this method. Otherwise, all bond orders will be
        set to an integer value in the building block.

        Parameters
        ----------
        molecule : :class:`rdkit.Mol`
            The molecule.

        functional_groups : :class:`iterable`, optional
            An :class:`iterable` of :class:`.FunctionalGroup` or
            :class:`.FunctionalGroupFactory` or both.
            :class:`.FunctionalGroup` instances are added to the
            building block and :class:`.FunctionalGroupFactory`
            instances are used to create :class:`.FunctionalGroup`
            instances the building block should hold.
            :class:`.FunctionalGroup` instances are used to identify
            which atoms are modified during
            :class:`.ConstructedMolecule` construction.

        placer_ids : :class:`tuple` of :class:`int`, optional
            The ids of *placer* atoms. These are the atoms which should
            be used for calculating the position of the building block.
            Depending on the values passed to `placer_ids`,
            and the functional groups in the building block, different
            *placer* ids will be used by the building block.

            #. `placer_ids` is passed to the initializer: the passed
               *placer* ids will be used by the building block.

            #. `placer_ids` is ``None`` and the building block has
               functional groups: The *placer* ids of the functional
               groups will be used as the *placer* ids of the building
               block.

            #. `placer_ids` is ``None`` and `functional_groups` is
               empty. All atoms of the molecule will be used for
               *placer* ids.

        Returns
        -------
        :class:`BuildingBlock`
            The molecule.

        """

        building_block = cls.__new__(cls)
        building_block._init_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
            placer_ids=placer_ids,
        )
        return building_block

    def _init_from_rdkit_mol(
        self,
        molecule,
        functional_groups,
        placer_ids,
    ):
        """
        Initialize from an :mod:`rdkit` molecule.

        Parameters
        ----------
        molecule : :class:`rdkit.Mol`
            The molecule.

        functional_groups : :class:`iterable`
            An :class:`iterable` of :class:`.FunctionalGroup` or
            :class:`.FunctionalGroupFactory` or both.
            :class:`.FunctionalGroup` instances are added to the
            building block and :class:`.FunctionalGroupFactory`
            instances are used to create :class:`.FunctionalGroup`
            instances the building block should hold.
            :class:`.FunctionalGroup` instances are used to identify
            which atoms are modified during
            :class:`.ConstructedMolecule` construction.

        placer_ids : :class:`tuple` of :class:`int`
            The ids of *placer* atoms. These are the atoms which should
            be used for calculating the position of the building block.
            Depending on the values passed to `placer_ids`,
            and the functional groups in the building block, different
            *placer* ids will be used by the building block.

            #. `placer_ids` is passed to the initializer: the passed
               *placer* ids will be used by the building block.

            #. `placer_ids` is ``None`` and the building block has
               functional groups: The *placer* ids of the functional
               groups will be used as the *placer* ids of the building
               block.

            #. `placer_ids` is ``None`` and `functional_groups` is
               empty. All atoms of the molecule will be used for
               *placer* ids.

        Returns
        -------
        None : :class:`NoneType`

        """

        atoms = tuple(
            Atom(a.GetIdx(), a.GetAtomicNum(), a.GetFormalCharge())
            for a in molecule.GetAtoms()
        )
        bonds = tuple(
            Bond(
                atom1=atoms[b.GetBeginAtomIdx()],
                atom2=atoms[b.GetEndAtomIdx()],
                order=(
                    9 if b.GetBondType() == rdkit.BondType.DATIVE
                    else b.GetBondTypeAsDouble()
                )
            )
            for b in molecule.GetBonds()
        )
        position_matrix = molecule.GetConformer().GetPositions()

        super().__init__(atoms, bonds, position_matrix)
        self._with_functional_groups(self._extract_functional_groups(
            functional_groups=functional_groups,
        ))
        self._placer_ids = self._normalize_placer_ids(
            placer_ids=placer_ids,
            functional_groups=self._functional_groups,
        )
        self._core_ids = frozenset(self._get_core_ids(
            functional_groups=self._functional_groups,
        ))

    def _normalize_placer_ids(
        self,
        placer_ids: Optional[tuple[int, ...]],
        functional_groups: Collection[FunctionalGroup],
    ) -> frozenset[int]:
        """
        Get the final *placer* ids.

        Parameters:

            placer_ids: The ids of *placer* atoms or ``None``.

            functional_groups:
                The :class:`.FunctionalGroup` instances of the building
                block.

        Returns:

            Depending on the input values, this function will return
            different things.

            #. `placer_ids` is a :class:`tuple` of :class`int`: the
                `placer_ids` will be returned.

            #. `placer_ids` is ``None`` and `functional_groups` is not
                empty: The *placer* ids of the functional groups will
                be returned.

            #. `placer_ids` is ``None`` and `functional_groups` is
               empty. The ids of all atoms in the building block will
               be returned.

        """

        if placer_ids is not None:
            return frozenset(placer_ids)

        if functional_groups:
            return frozenset(flatten(
                functional_group.get_placer_ids()
                for functional_group in functional_groups
            ))

        return frozenset(atom.get_id() for atom in self._atoms)

    def _get_core_ids(
        self,
        functional_groups: Iterable[FunctionalGroup],
    ) -> Iterator[int]:
        """
        Get the final *core* ids.

        This method may return duplicate ids.

        Parameters:

            functional_groups:
                The :class:`.FunctionalGroup` instances of the building
                block.

        Yields:

            The id of an atom defining the core of the molecule.

        """

        functional_group_atom_ids = {
            atom_id
            for functional_group in functional_groups
            for atom_id in functional_group.get_atom_ids()
        }
        for atom in self._atoms:
            atom_id = atom.get_id()
            if atom_id not in functional_group_atom_ids:
                yield atom_id

        for functional_group in functional_groups:
            for atom_id in functional_group.get_core_atom_ids():
                yield atom_id

    def _extract_functional_groups(self, functional_groups):
        """
        Yield functional groups.

        The input can be a mixture of :class:`.FunctionalGroup` and
        :class:`.FunctionalGroupFactory`. The output yields
        :class:`.FunctionalGroup` instances only. Either those
        held directly in `functional_groups` or created by the
        factories in `functional_groups`.

        Parameters
        ----------
        functional_groups : :class:`iterable`
            Can be an :class:`iterable` of both
            :class:`.FunctionalGroup` and
            :class:`.FunctionalGroupFactory`.

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group from `functional_groups`, or created
            by a factory in `functional_groups`.

        """

        for functional_group in functional_groups:
            if isinstance(functional_group, FunctionalGroup):
                yield functional_group
            else:
                # Else it's a factory.
                yield from functional_group.get_functional_groups(self)

    def _with_functional_groups(self, functional_groups):
        """
        Modify the molecule.

        """

        self._functional_groups = tuple(functional_groups)
        return self

    def with_functional_groups(self, functional_groups):
        """
        Return a clone with specific functional groups.

        Parameters
        ----------
        functional_groups : :class:`iterable`
            :class:`.FunctionalGroup` instances which the clone
            should have.

        Returns
        -------
        :class:`.BuildingBlock`
            The clone. Has the same type as the original molecule.

        """

        return self.clone()._with_functional_groups(functional_groups)

    def _with_canonical_atom_ordering(self):
        ordering = rdkit.CanonicalRankAtoms(self.to_rdkit_mol())
        super()._with_canonical_atom_ordering()
        id_map = {
            old_id: new_id
            for old_id, new_id in enumerate(ordering)
        }
        self._functional_groups = tuple(
            functional_group.with_ids(id_map)
            for functional_group in self._functional_groups
        )
        self._placer_ids = frozenset(
            id_map[placer_id]
            for placer_id in self._placer_ids
        )
        self._core_ids = frozenset(
            id_map[core_id]
            for core_id in self._core_ids
        )
        return self

    def get_num_functional_groups(self):
        """
        Return the number of functional groups.

        Returns
        -------
        :class:`int`
            The number of functional groups in the building block.

        """

        return len(self._functional_groups)

    def get_functional_groups(
        self,
        fg_ids: Optional[Union[int, Iterable[int]]] = None,
    ) -> Iterable[FunctionalGroup]:
        """
        Yield the functional groups, ordered by id.

        Parameters:

            fg_ids:
                The ids of functional groups yielded. If ``None``, then
                all functional groups are yielded. Can be a single
                :class:`int`, if a single functional group is
                desired.

        Yields:

            A functional group of the building block.

        """

        if fg_ids is None:
            fg_ids = range(len(self._functional_groups))
        elif isinstance(fg_ids, int):
            fg_ids = (fg_ids, )

        for fg_id in fg_ids:
            yield self._functional_groups[fg_id]

    def clone(self):
        clone = super().clone()
        clone._functional_groups = self._functional_groups
        clone._placer_ids = self._placer_ids
        clone._core_ids = self._core_ids
        return clone

    def get_num_placers(self) -> int:
        """
        Return the number of *placer* atoms in the building block.

        Returns:

            The number of *placer* atoms in the building block.

        """

        return len(self._placer_ids)

    def get_placer_ids(self):
        """
        Yield the ids of *placer* atoms.

        *Placer* atoms are those, which should be used to calculate
        the position of the building block.

        Yields
        ------
        :class:`int`
            The id of a *placer* atom.

        See Also
        --------
        :meth:`.FunctionalGroup.get_placer_ids`

        """

        yield from self._placer_ids

    def get_core_atom_ids(self):
        """
        Yield ids of atoms which form the core of the building block.

        This includes all atoms in the building block not part of a
        functional group, as well as any atoms in a functional group,
        specifically labelled as core atoms.

        Yields
        ------
        :class:`int`
            The id of a *core* atom.

        See Also
        --------
        :meth:`.FunctionalGroup.get_core_atom_ids`

        """

        yield from self._core_ids

    def __str__(self):
        if self._functional_groups:
            fg_repr = f', {self._functional_groups!r}'
        else:
            fg_repr = ''

        smiles = rdkit.MolToSmiles(
            mol=rdkit.RemoveHs(self.to_rdkit_mol()),
        )
        return f'{self.__class__.__name__}({smiles!r}{fg_repr})'

    def __repr__(self):
        return str(self)
