from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import RingAmine


class RingAmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.RingAmine` functional groups.

    """

    _functional_group_smarts = '[N]([H])([H])[#6]~[#6]([H])~[#6R1]'
    _bonder_smarts = [
        '[$([N]([H])([H])[#6]~[#6]([H])~[#6R1])]',
        '[$([#6]([H])(~[#6R1])~[#6][N]([H])[H])]',
    ]
    _deleter_smarts = (
            ['[$([H][N]([H])[#6]~[#6]([H])~[#6R1])]']*2 +
            ['[$([H][#6](~[#6R1])~[#6][N]([H])[H])]']
    )

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield RingAmine(
                atom1=atoms[1],
                bromine1=atoms[0],
                atom2=atoms[2],
                bromine2=atoms[3],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
