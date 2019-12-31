from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Iodo


class IodoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Iodo` instances.

    """

    _functional_group_smarts = '*[I]'
    _bonder_smarts = ['[$(*[I])]']
    _deleter_smarts = ['[$([I]*)]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Iodo(
                iodine=atoms[1],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )
