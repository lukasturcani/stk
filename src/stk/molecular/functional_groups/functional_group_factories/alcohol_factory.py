from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Alcohol


class AlcoholFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Alcohol` instances.

    """

    _functional_group_smarts = '[*][O][H]'
    _bonder_smarts = ['[$([O][H])]']
    _deleter_smarts = ['[$([H][O])]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Alcohol(
                oxygen=atoms[1],
                hydrogen=atoms[2],
                atom=atoms[0],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
