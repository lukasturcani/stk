from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Thioacid


class ThioacidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Thioacid` instances.

    """

    _functional_group_smarts = '[*][C](=[O])[S][H]'
    _bonder_smarts = ['[$([C](=[O])[S][H])]']
    _deleter_smarts = [
        '[$([H][S][C](=[O]))]',
        '[$([S]([H])[C](=[O]))]',
    ]

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Thioacid(
                carbon=atoms[1],
                oxygen=atoms[2],
                sulfur=atoms[3],
                hydrogen=atoms[4],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )
