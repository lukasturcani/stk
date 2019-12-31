from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Thiol


class ThiolFactory(SmartsFunctionalGroupFactory):
    _functional_group_smarts = '[*][S][H]'
    _bonder_smarts = ['[$([S][H])]']
    _deleter_smarts = ['[$([H][S])]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Thiol(
                sulfur=atoms[1],
                hydrogen=atoms[2],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )
