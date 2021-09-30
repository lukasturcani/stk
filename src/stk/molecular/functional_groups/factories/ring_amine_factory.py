"""
Ring Amine Factory
==================

"""

from __future__ import annotations

import typing
from collections import abc

from .functional_group_factory import FunctionalGroupFactory
from .utilities import get_atom_ids
from ..functional_groups import RingAmine
from ...molecule import Molecule
from ...atoms import elements


class RingAmineFactory(FunctionalGroupFactory):
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
        molecule: Molecule,
    ) -> abc.Iterable[RingAmine]:

        ids = get_atom_ids(
            query='[N]([H])([H])[#6]~[#6]([H])~[#6R1]',
            molecule=molecule,
        )
        for atom_ids in ids:
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield RingAmine(
                nitrogen=typing.cast(elements.N, atoms[0]),
                hydrogen1=typing.cast(elements.H, atoms[1]),
                hydrogen2=typing.cast(elements.H, atoms[2]),
                carbon1=typing.cast(elements.C, atoms[3]),
                carbon2=typing.cast(elements.C, atoms[4]),
                hydrogen3=typing.cast(elements.H, atoms[5]),
                carbon3=typing.cast(elements.C, atoms[6]),
            )
