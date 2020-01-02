from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ...functional_groups import Difluoro


class DifluoroFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Difluoro` instances.

    """

    def __init__(self, bonders=(1, 2), deleters=(0, 3)):
        super().__init__('[F][#6]~[#6][F]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Difluoro(
                atom1=atoms[1],
                fluorine1=atoms[0],
                atom2=atoms[2],
                fluorine2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
