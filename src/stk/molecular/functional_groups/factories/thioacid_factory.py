from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Thioacid


class ThioacidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Thioacid` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(3, 4)):
        super().__init__('[*][C](=[O])[S][H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Thioacid(
                carbon=atoms[1],
                oxygen=atoms[2],
                sulfur=atoms[3],
                hydrogen=atoms[4],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
