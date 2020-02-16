from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import CarboxylicAcid


class CarboxylicAcidFactory(FunctionalGroupFactory):
    """
    Creates :class:`.CarboxylicAcid` instances.

    Creates functional groups from substructures, which match the
    ``[*][C](=[O])[O][H]`` functional group string.

    """

    def __init__(self, bonders=(1, ), deleters=(3, 4)):
        """
        Initialize a :class:`.CarboxylicAcidFactory` instance.

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
        for atom_ids in _get_atom_ids('[*][C](=[O])[O][H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield CarboxylicAcid(
                carbon=atoms[1],
                oxygen1=atoms[2],
                oxygen2=atoms[3],
                hydrogen=atoms[4],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
