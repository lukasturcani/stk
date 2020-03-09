"""
Aldehyde Factory
================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Aldehyde


class AldehydeFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Aldehyde` instances.

    Creates functional groups from substructures, which match the
    ``[*][C](=[O])[H]`` functional group string.

    """

    def __init__(self, bonders=(1, ), deleters=(2, )):
        """
        Initialize a :class:`.AldehydeFactory` instance.

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
        for atom_ids in _get_atom_ids('[*][C](=[O])[H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Aldehyde(
                carbon=atoms[1],
                oxygen=atoms[2],
                hydrogen=atoms[3],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
