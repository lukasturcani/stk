from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import BoronicAcid


class BoronicAcidFactory(FunctionalGroupFactory):
    """
    Creates :class:`.BoronicAcid` instances.

    Creates functional groups from substructures, which match the
    ``[*][B]([O][H])[O][H]`` functional group string.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3, 4, 5)):
        """
        Initialize a :class:`.BoronicAcidFactory` instance.

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
        ids = _get_atom_ids('[*][B]([O][H])[O][H]', molecule)
        for atom_ids in ids:
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield BoronicAcid(
                boron=atoms[1],
                oxygen1=atoms[2],
                hydrogen1=atoms[3],
                oxygen2=atoms[4],
                hydrogen2=atoms[5],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
