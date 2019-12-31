from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Amide


class AmideFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Amide` instances.

    """

    _functional_group_smarts = '[*][C](=[O])[N]([H])[H]'
    _bonder_smarts = ['[$([C](=[O])[N]([H])[H])]']
    _deleter_smarts = (
        ['[$([N]([H])([H])[C](=[O]))]'] +
        ['[$([H][N]([H])[C](=[O]))]']*2
    )

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Amide(
                carbon=atoms[1],
                oxygen=atoms[2],
                nitrogen=atoms[3],
                hydrogen1=atoms[4],
                hydrogen2=atoms[5],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )
