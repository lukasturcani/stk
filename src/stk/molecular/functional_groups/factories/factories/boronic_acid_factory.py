from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import BoronicAcid


class BoronicAcidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.BoronicAcid` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3, 4, 5)):
        super().__init__('[*][B]([O][H])[O][H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
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
