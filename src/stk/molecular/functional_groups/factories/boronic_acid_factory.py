"""
Boronic Acid Factory
====================

"""

from __future__ import annotations

import typing
from collections import abc

from .functional_group_factory import FunctionalGroupFactory
from .utilities import get_atom_ids
from ..functional_groups import BoronicAcid
from ...molecule import Molecule
from ...elements import O, H, B


__all__ = (
    'BoronicAcidFactory',
)


_ValidIndex = typing.Literal[0, 1, 2, 3, 4, 5]


class BoronicAcidFactory(FunctionalGroupFactory):
    """
    Creates :class:`.BoronicAcid` instances.

    Creates functional groups from substructures, which match the
    ``[*][B]([O][H])[O][H]`` functional group string.

    Examples:


        *Creating Functional Groups with the Factory*

        You want to create a building block which has
        :class:`.BoronicAcid` functional groups. You want the boron
        atom in those functional groups to be the *bonder* atom and the
        OH groups to be *deleter* atoms.

        .. testcode:: creating-functional-groups-with-the-factory

            import stk

            building_block = stk.BuildingBlock(
                smiles='OB(O)CCCB(O)O',
                functional_groups=(stk.BoronicAcidFactory(), ),
            )

        .. testcode:: creating-functional-groups-with-the-factory
            :hide:

            assert all(
                isinstance(functional_group, stk.BoronicAcid)
                for functional_group
                in building_block.get_functional_groups()
            )
            assert building_block.get_num_functional_groups() == 2

        *Changing the Bonder and Deleter Atoms*

        You want to create a building block which has
        :class:`.BoronicAcid` functional groups. You want the oxygen
        atoms to be treated as *bonder* atoms, and the hydrogen atoms
        to be treated as *deleter* atoms.

        .. testcode:: changing-the-bonder-and-deleter-atoms

            import stk

            boronic_acid_factory = stk.BoronicAcidFactory(
                # The indices of the oxygen atoms in the functional
                # group string (see docstring) are 2 and 4.
                bonders=(2, 4),
                # The indices of the hydrogen atoms in the
                # functional group string (see docstring) are 3 and 5.
                deleters=(3, 5),
            )
            building_block = stk.BuildingBlock(
                smiles='OB(O)CCCB(O)O',
                functional_groups=(boronic_acid_factory, ),
            )

        .. testcode:: changing-the-bonder-and-deleter-atoms
            :hide:

            fg1, fg2 = building_block.get_functional_groups()
            assert fg1.get_num_bonders() == 2
            assert sum(1 for _ in fg1.get_deleters()) == 2
            assert fg2.get_num_bonders() == 2
            assert sum(1 for _ in fg2.get_deleters()) == 2

            assert all(
                isinstance(atom, stk.O)
                for functional_group
                in building_block.get_functional_groups()
                for atom
                in functional_group.get_bonders()
            )
            assert all(
                isinstance(atom, stk.H)
                for functional_group
                in building_block.get_functional_groups()
                for atom
                in functional_group.get_deleters()
            )

    See Also:

        :class:`.GenericFunctionalGroup`
            Defines *bonders* and  *deleters*.

    """

    def __init__(
        self,
        bonders: tuple[_ValidIndex, ...] = (1, ),
        deleters: tuple[_ValidIndex, ...] = (2, 3, 4, 5),
        placers: typing.Optional[tuple[_ValidIndex, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.BoronicAcidFactory` instance.

        Parameters:

            bonders:
                The indices of atoms in the functional group string,
                which are *bonder* atoms.

            deleters:
                The indices of atoms in the functional group string,
                which are *deleter* atoms.

            placers:
                The indices of atoms in the functional group string,
                which are *placer* atoms. If ``None``, `bonders` will
                be used.

        """

        self._bonders = bonders
        self._deleters = deleters
        self._placers = bonders if placers is None else placers

    def get_functional_groups(
        self,
        molecule: Molecule,
    ) -> abc.Iterable[BoronicAcid]:

        for atom_ids in get_atom_ids(
            query='[*][B]([O][H])[O][H]',
            molecule=molecule,
        ):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield BoronicAcid(
                boron=typing.cast(B, atoms[1]),
                oxygen1=typing.cast(O, atoms[2]),
                hydrogen1=typing.cast(H, atoms[3]),
                oxygen2=typing.cast(O, atoms[4]),
                hydrogen2=typing.cast(H, atoms[5]),
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
