from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ...functional_groups import Amide


class AmideFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Amide` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(3, 4, 5)):
        super().__init__('[*][C](=[O])[N]([H])[H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Amide(
                carbon=atoms[1],
                oxygen=atoms[2],
                nitrogen=atoms[3],
                hydrogen1=atoms[4],
                hydrogen2=atoms[5],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
