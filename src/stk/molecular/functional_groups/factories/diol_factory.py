"""
Diol Factory
============

"""

from __future__ import annotations

import typing
from collections import abc

from . import functional_group_factory as _functional_group_factory
from . import utilities as _utilities
from .. import functional_groups as _functional_groups
from ... import molecule as _molecule
from ...atoms import elements as _elements

__all__ = (
    'DiolFactory',
)

_ValidIndex = typing.Literal[0, 1, 2, 3, 4, 5]


class DiolFactory(
    _functional_group_factory.FunctionalGroupFactory,
):
    """
    Creates :class:`.Diol` instances.

    Creates functional groups from substructures, which match the
    ``[H][O][#6]~[#6][O][H]`` functional group string.

    Examples:

        *Creating Functional Groups with the Factory*

        You want to create a building block which has :class:`.Diol`
        functional groups. You want the carbon atoms in those
        functional groups to be the *bonder* atoms, and the OH groups
        to be a leaving groups.

        .. testcode:: creating-functional-groups-with-the-factory

            import stk

            building_block = stk.BuildingBlock(
                smiles='CCCC(O)C(O)CCCC',
                functional_groups=(stk.DiolFactory(), ),
            )

        .. testcode:: creating-functional-groups-with-the-factory
            :hide:

            assert all(
                isinstance(functional_group, stk.Diol)
                for functional_group
                in building_block.get_functional_groups()
            )
            assert building_block.get_num_functional_groups() == 1

        *Changing the Bonder and Deleter Atoms*

        You want to create a building block which has :class:`.Diol`
        functional groups. You want the oxygen atoms to be the *bonder*
        atoms and the hydrogen atoms to be the *deleter* atoms.

        .. testcode:: changing-the-bonder-and-deleter-atoms

            import stk

            diol_factory = stk.DiolFactory(
                # The indices of the oxygen atoms in the functional
                # group string (see docstring) are 1 and 4.
                bonders=(1, 4),
                # The indices of the hydrogen atoms in the functional
                # group string (see docstring) are 0 and 5.
                deleters=(0, 5),
            )
            building_block = stk.BuildingBlock(
                smiles='CCCC(O)C(O)CCCC',
                functional_groups=(diol_factory, ),
            )

        .. testcode:: changing-the-bonder-and-deleter-atoms
            :hide:

            fg, = building_block.get_functional_groups()
            assert fg.get_num_bonders() == 2
            assert sum(1 for _ in fg.get_deleters()) == 2

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
        bonders: tuple[_ValidIndex, ...] = (2, 3),
        deleters: tuple[_ValidIndex, ...] = (0, 1, 4, 5),
        placers: typing.Optional[tuple[_ValidIndex, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.DiolFactory` instance.

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
        molecule: _molecule.Molecule,
    ) -> abc.Iterable[_functional_groups.Diol]:

        for atom_ids in _utilities.get_atom_ids(
            query='[H][O][#6]~[#6][O][H]',
            molecule=molecule,
        ):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield _functional_groups.Diol(
                hydrogen1=typing.cast(_elements.H, atoms[0]),
                oxygen1=typing.cast(_elements.O, atoms[1]),
                atom1=atoms[2],
                atom2=atoms[3],
                oxygen2=typing.cast(_elements.O, atoms[4]),
                hydrogen2=typing.cast(_elements.H, atoms[5]),
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
