from .smarts_functioanal_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import CarboxylicAcid


class CarboxylicAcidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.CarboxylicAcid` instances.

    """

    _functional_group_smarts = '[*][C](=[O])[O][H]'
    _bonder_smarts = ['[$([C](=[O])[O][H])]']
    _deleter_smarts = [
        '[$([H][O][C](=[O]))]',
        '[$([O]([H])[C](=[O]))]',
    ]

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield CarboxylicAcid(
                carbon=atoms[1],
                oxygen1=atoms[2],
                oxygen2=atoms[3],
                hydrogen=atoms[4],
                atom=atoms[0],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
