from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Dibromo


class DibromoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Dibromo` instances.

    """

    def __init__(self, bonders=(1, 2), deleters=(0, 3)):
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
