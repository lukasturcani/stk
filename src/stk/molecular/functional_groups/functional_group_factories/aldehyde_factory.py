from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Aldehyde


class AldehydeFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Aldehyde` instances.

    """

    _functional_group_smarts = '[*][C](=[O])[H]'
    _bonder_smarts = ['[$([C](=[O])[H])]']
    _deleter_smarts = ['[$([O]=[C][H])]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Aldehyde(
                carbon=atoms[1],
                oxygen=atoms[2],
                hydrogen=atoms[3],
                atom=atoms[0],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
