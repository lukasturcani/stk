"""
Aldehyde Factory
================

"""

from __future__ import annotations

import typing
from collections import abc

from .functional_group_factory import FunctionalGroupFactory
from .utilities import get_atom_ids
from ..functional_groups import Aldehyde
from ...molecule import Molecule
from ...atoms import C, O, H


__all__ = (
    'AldehydeFactory',
)


_ValidIndex = typing.Literal[0, 1, 2, 3]


class AldehydeFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Aldehyde` instances.

    Creates functional groups from substructures, which match the
    ``[*][C](=[O])[H]`` functional group string.

    Examples:

        *Creating Functional Groups with the Factory*

        You want to create a building block which has
        :class:`.Aldehyde` functional groups. You want the carbon atom
        in those functional groups to be the bonder atom, and the
        oxygen atom to be the deleter atom.

        .. testcode:: creating-functional-groups-with-the-factory

            import stk

            building_block = stk.BuildingBlock(
                smiles='O=CCC=O',
                functional_groups=(stk.AldehydeFactory(), ),
            )

        .. testcode:: creating-functional-groups-with-the-factory
            :hide:

            assert all(
                isinstance(functional_group, stk.Aldehyde)
                for functional_group
                in building_block.get_functional_groups()
            )
            assert building_block.get_num_functional_groups() == 2

        *Changing the Bonder and Deleter Atoms*

        You want to create a building block which has
        :class:`.Aldehyde` functional groups. You want the carbon atom
        to be the bonder atom and the hydrogen atom to be the deleter
        atom.

        .. testcode:: changing-the-bonder-and-deleter-atoms

            import stk

            aldehyde_factory = stk.AldehydeFactory(
                # The index of the carbon atom in the functional
                # group string (see docstring) is 1.
                bonders=(1, ),
                # The index of the hydrogen atom in the functional
                # group string (see docstring) is 3.
                deleters=(3, ),
            )
            building_block = stk.BuildingBlock(
                smiles='O=CCC=O',
                functional_groups=(aldehyde_factory, ),
            )

        .. testcode:: changing-the-bonder-and-deleter-atoms
            :hide:

            fg1, fg2 = building_block.get_functional_groups()
            assert fg1.get_num_bonders() == 1
            assert sum(1 for _ in fg1.get_deleters()) == 1
            assert fg2.get_num_bonders() == 1
            assert sum(1 for _ in fg2.get_deleters()) == 1

            assert all(
                isinstance(atom, stk.C)
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
        deleters: tuple[_ValidIndex, ...] = (2, ),
        placers: typing.Optional[tuple[_ValidIndex, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.AldehydeFactory` instance.

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
    ) -> abc.Iterable[Aldehyde]:

        for atom_ids in get_atom_ids(
            query='[*][C](=[O])[H]',
            molecule=molecule,
        ):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Aldehyde(
                carbon=typing.cast(C, atoms[1]),
                oxygen=typing.cast(O, atoms[2]),
                hydrogen=typing.cast(H, atoms[3]),
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
