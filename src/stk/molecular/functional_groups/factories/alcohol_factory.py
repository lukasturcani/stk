"""
Alcohol Factory
===============

"""

from __future__ import annotations

import typing
from collections import abc

from .functional_group_factory import FunctionalGroupFactory
from .utilities import get_atom_ids
from ..functional_groups import Alcohol
from ...molecule import Molecule
from ...atoms import O, H


__all__ = (
    'AlcoholFactory',
)


_ValidIndex = typing.Literal[0, 1, 2]


class AlcoholFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Alcohol` instances.

    Creates functional groups from substructures, which match the
    ``[*][O][H]`` functional group string.

    Examples:

        *Creating Functional Groups with the Factory*

        You want to create a building block which has :class:`.Alcohol`
        functional groups. You want the oxygen atom in those functional
        groups to be the bonder atom, and the hydrogen atom to be the
        deleter atom.

        .. testcode:: creating-functional-groups-with-the-factory

            import stk

            building_block = stk.BuildingBlock(
                smiles='OCCCO',
                functional_groups=(stk.AlcoholFactory(), ),
            )

        .. testcode:: creating-functional-groups-with-the-factory
            :hide:

            assert all(
                isinstance(functional_group, stk.Alcohol)
                for functional_group
                in building_block.get_functional_groups()
            )
            assert building_block.get_num_functional_groups() == 2

        *Changing the Bonder and Deleter Atoms*

        You want to create a building block which has :class:`.Alcohol`
        functional groups. You want the OH group to be
        treated as a leaving group. This means the non-hydrogen atom
        bonded to oxygen is the *bonder* atom and both the oxygen and
        hydrogen atoms are *deleter* atoms.

        .. testcode:: changing-the-bonder-and-deleter-atoms

            import stk

            alcohol_factory = stk.AlcoholFactory(
                # The index of the non-hydrogen atom connected to
                # oxygen is 0 in the functional group string (see
                # docstring).
                bonders=(0, ),
                # The indices of the oxygen and hydrogen atoms in the
                # functional group string (see docstring) are
                # 1 and 2, respectively.
                deleters=(1, 2),
            )
            building_block = stk.BuildingBlock(
                smiles='OCCCO',
                functional_groups=(alcohol_factory, ),
            )

        .. testcode:: changing-the-bonder-and-deleter-atoms
            :hide:

            fg1, fg2 = building_block.get_functional_groups()
            assert fg1.get_num_bonders() == 1
            assert sum(1 for _ in fg1.get_deleters()) == 2
            assert fg2.get_num_bonders() == 1
            assert sum(1 for _ in fg2.get_deleters()) == 2

            assert all(
                isinstance(atom, stk.C)
                for functional_group
                in building_block.get_functional_groups()
                for atom
                in functional_group.get_bonders()
            )
            assert all(
                isinstance(atom, (stk.O, stk.H))
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
        deleters: tuple[_ValidIndex, ...] = (2, ),
        placers: typing.Optional[tuple[_ValidIndex, ...]] = None,
    ) -> None:
        """
        Initialize an :class:`.AlcoholFactory` instance.

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
    ) -> abc.Iterable[Alcohol]:

        for atom_ids in get_atom_ids('[*][O][H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Alcohol(
                oxygen=typing.cast(O, atoms[1]),
                hydrogen=typing.cast(H, atoms[2]),
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
