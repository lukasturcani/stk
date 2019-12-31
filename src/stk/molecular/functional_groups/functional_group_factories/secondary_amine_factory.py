from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import SecondaryAmine


class SecondaryAmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.SecondaryAmine` instances.

    """

    _functional_group_smarts = '[H][N]([#6])[#6]'
    _bonder_smarts = ['[$([N]([H])([#6])[#6])]']
    _deleter_smarts = ['[$([H][N]([#6])[#6])]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield SecondaryAmine(
                nitrogen=atoms[1],
                hydrogen=atoms[0],
                atom1=atoms[2],
                atom2=atoms[3],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
