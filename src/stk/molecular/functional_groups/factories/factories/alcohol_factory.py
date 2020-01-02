from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Alcohol


class AlcoholFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Alcohol` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(2, )):
        super().__init__('[*][O][H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Alcohol(
                oxygen=atoms[1],
                hydrogen=atoms[2],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
