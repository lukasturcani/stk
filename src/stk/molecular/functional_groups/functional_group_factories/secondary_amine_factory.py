from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import SecondaryAmine


class SecondaryAmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.SecondaryAmine` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(0, )):
        super().__init__('[H][N]([#6])[#6]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield SecondaryAmine(
                nitrogen=atoms[1],
                hydrogen=atoms[0],
                atom1=atoms[2],
                atom2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
