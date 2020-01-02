from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ...functional_groups import CarboxylicAcid


class CarboxylicAcidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.CarboxylicAcid` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(3, 4)):
        super().__init__('[*][C](=[O])[O][H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
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
