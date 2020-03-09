"""
Secondary Amino Factory
=======================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import SecondaryAmino


class SecondaryAminoFactory(FunctionalGroupFactory):
    """
    Creates :class:`.SecondaryAmino` instances.

    Creates functional groups from substructures, which match the
    ``[H][N]([#6])[#6]`` functional group string.

    """

    def __init__(self, bonders=(1, ), deleters=(0, )):
        """
        Initialize a :class:`.SecondaryAminoFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are bonder atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are deleter atoms.

        """

        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids('[H][N]([#6])[#6]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield SecondaryAmino(
                nitrogen=atoms[1],
                hydrogen=atoms[0],
                atom1=atoms[2],
                atom2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
