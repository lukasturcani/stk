"""
Ring Amine Factory
==================

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
    'RingAmineFactory',
)


class RingAmineFactory(
    _functional_group_factory.FunctionalGroupFactory,
):
    """
    Creates :class:`.RingAmine` functional groups.

    Creates functional groups from substructures, which match the
    ``[N]([H])([H])[#6]~[#6]([H])~[#6R1]`` functional group string.

    """

    # Keep this __init__() method, which does nothing, here for the
    # docstring.
    def __init__(self):
        """
        Initialize a :class:`RingAmineFactory` instance.

        """

        return super().__init__()

    def get_functional_groups(
        self,
        molecule: _molecule.Molecule,
    ) -> abc.Iterable[_functional_groups.RingAmine]:

        for atom_ids in _utilities.get_atom_ids(
            query='[N]([H])([H])[#6]~[#6]([H])~[#6R1]',
            molecule=molecule,
        ):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield _functional_groups.RingAmine(
                nitrogen=typing.cast(_elements.N, atoms[0]),
                hydrogen1=typing.cast(_elements.H, atoms[1]),
                hydrogen2=typing.cast(_elements.H, atoms[2]),
                carbon1=typing.cast(_elements.C, atoms[3]),
                carbon2=typing.cast(_elements.C, atoms[4]),
                hydrogen3=typing.cast(_elements.H, atoms[5]),
                carbon3=typing.cast(_elements.C, atoms[6]),
            )
