from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Alkyne


class TerminalAlkyneFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Alkyne` instances.

    Creates functional groups from substructures, which match the
    ``[*][C]#[C][H]`` functional group string.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3)):
        """
        Initialize a :class:`.TerminalAlkyneFactory` instance.

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
        for atom_ids in _get_atom_ids('[*][C]#[C][H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Alkyne(
                atom1=atoms[0],
                carbon1=atoms[1],
                carbon2=atoms[2],
                atom2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
