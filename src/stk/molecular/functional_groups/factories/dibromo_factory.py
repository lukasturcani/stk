from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Dibromo


class DibromoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Dibromo` instances.

    Creates functional groups from substructures, which match the
    ``[Br][#6]~[#6][Br]`` functional group string.

    """

    def __init__(self, bonders=(1, 2), deleters=(0, 3)):
        """
        Initialize a :class:`.DibromoFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are bonder atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are deleter atoms.

        """

        super().__init__('[Br][#6]~[#6][Br]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Dibromo(
                atom1=atoms[1],
                bromine1=atoms[0],
                atom2=atoms[2],
                bromine2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
