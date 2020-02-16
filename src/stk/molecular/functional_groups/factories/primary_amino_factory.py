from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import PrimaryAmino


class PrimaryAminoFactory(FunctionalGroupFactory):
    """
    Creates :class:`.PrimaryAmino` instances.

    Creates functional groups from substructures, which match the
    ``[*][N]([H])[H]`` functional group string.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3)):
        """
        Initialize a :class:`.PrimaryAminoFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are bonder atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are deleter atoms.

        """

        """
        Initialize an :class:`.AmineFactory`.

        """

        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids('[*][N]([H])[H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield PrimaryAmino(
                nitrogen=atoms[1],
                hydrogen1=atoms[2],
                hydrogen2=atoms[3],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
