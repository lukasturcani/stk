"""
SMARTS Functional Group Factory
===============================

"""

from __future__ import annotations

from typing import Optional
from collections import abc

from . import functional_group_factory as _functional_group_factory
from . import utilities as _utilities
from .. import functional_groups as _functional_groups
from ... import molecule as _molecule


__all__ = (
    'SmartsFunctionalGroupFactory',
)


class SmartsFunctionalGroupFactory(
    _functional_group_factory.FunctionalGroupFactory,
):
    """
    Creates :class:`.GenericFunctionalGroup` instances.

    Examples:

        *Using SMARTS to Define Functional Groups*

        You want to create a building block which has
        :class:`.GenericFunctionalGroup` functional groups based on the
        SMARTS string: ``[Br][C]``.
        You want the ``C`` atom to be the *bonder* atom, and the
        ``Br`` atom to be the *deleter* atom.

        .. testcode:: using-smarts-to-define-functional-groups

            import stk

            building_block = stk.BuildingBlock(
                smiles='BrCCCBr',
                functional_groups=(
                    stk.SmartsFunctionalGroupFactory(
                        smarts='[Br][C]',
                        bonders=(1, ),
                        deleters=(0, ),
                    ),
                ),
            )

        .. testcode:: using-smarts-to-define-functional-groups
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
                isinstance(atom, stk.Br)
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
        smarts: str,
        bonders: tuple[int, ...],
        deleters: tuple[int, ...],
        placers: Optional[tuple[int, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.SmartsFunctionalGroupFactory` instance.

        Parameters:

            smarts:
                The SMARTS defining the functional group.

            bonders:
                The indices of atoms in `smarts`, which are *bonder*
                atoms.

            deleters:
                The indices of atoms in `smarts`, which are *deleter*
                atoms.

            placers:
                The indices of atoms in `smarts`, which are *placer*
                atoms. If ``None``, the *bonder* atoms will be used.

        """

        self._smarts = smarts
        self._bonders = bonders
        self._deleters = deleters
        self._placers = bonders if placers is None else placers

    def get_functional_groups(
        self,
        molecule: _molecule.Molecule,
    ) -> abc.Iterable[_functional_groups.GenericFunctionalGroup]:

        for atom_ids in _utilities.get_atom_ids(
            query=self._smarts,
            molecule=molecule,
        ):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield _functional_groups.GenericFunctionalGroup(
                atoms=atoms,
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
