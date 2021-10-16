"""
Amide Factory
=============

"""

from __future__ import annotations

import typing
from collections import abc

from .functional_group_factory import FunctionalGroupFactory
from .utilities import get_atom_ids
from ..functional_groups import Amide
from ...molecule import Molecule
from ...elements import C, O, H, N


__all__ = (
    'AmideFactory',
)


_ValidIndex = typing.Literal[0, 1, 2, 3, 4, 5]


class AmideFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Amide` instances.

    Creates functional groups from substructures, which match the
    ``[*][C](=[O])[N]([H])[H]`` functional group string.

    Examples:

        *Creating Functional Groups with the Factory*

        You want to create a building block which has :class:`.Amide`
        functional groups. You want the carbon atom in those functional
        groups to be the bonder atom, and the amino group to be a
        leaving group.

        .. testcode:: creating-functional-groups-with-the-factory

            import stk

            building_block = stk.BuildingBlock(
                smiles='NC(=O)CC(=O)N',
                functional_groups=(stk.AmideFactory(), ),
            )

        .. testcode:: creating-functional-groups-with-the-factory
            :hide:

            assert all(
                isinstance(functional_group, stk.Amide)
                for functional_group
                in building_block.get_functional_groups()
            )
            assert building_block.get_num_functional_groups() == 2

        *Changing the Bonder and Deleter Atoms*

        You want to create a building block which has :class:`.Amide`
        functional groups. You want the carbon atom to be the bonder
        atom and the oxygen atom to be the deleter atom.

        .. testcode:: changing-the-bonder-and-deleter-atoms

            import stk

            amide_factory = stk.AmideFactory(
                # The index of the carbon atom in the functional
                # group string (see docstring) is 1.
                bonders=(1, ),
                # The index of the oxygen atom in the functional
                # group string (see docstring) is 2.
                deleters=(2, ),
            )
            building_block = stk.BuildingBlock(
                smiles='NC(=O)CC(=O)N',
                functional_groups=(amide_factory, ),
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
                isinstance(atom, stk.O)
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
        deleters: tuple[_ValidIndex, ...] = (3, 4, 5),
        placers: typing.Optional[tuple[_ValidIndex, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.AmideFactory` instance.

        Parameters:

            bonders:
                The indices of atoms in the functional group string,
                which are bonder atoms.

            deleters:
                The indices of atoms in the functional group string,
                which are deleter atoms.

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
    ) -> abc.Iterable[Amide]:

        for atom_ids in get_atom_ids(
            query='[*][C](=[O])[N]([H])[H]',
            molecule=molecule,
        ):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Amide(
                carbon=typing.cast(C, atoms[1]),
                oxygen=typing.cast(O, atoms[2]),
                nitrogen=typing.cast(N, atoms[3]),
                hydrogen1=typing.cast(H, atoms[4]),
                hydrogen2=typing.cast(H, atoms[5]),
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
