from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Diol


class DiolFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Diol` instances.

    """

    def __init__(self, bonders=(2, 3), deleters=(0, 1, 4, 5)):
        super().__init__('[H][O][#6]~[#6][O][H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Diol(
                hydrogen1=atoms[0],
                oxygen1=atoms[1],
                atom1=atoms[2],
                atom2=atoms[3],
                oxygen2=atoms[4],
                hydrogen2=atoms[5],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
