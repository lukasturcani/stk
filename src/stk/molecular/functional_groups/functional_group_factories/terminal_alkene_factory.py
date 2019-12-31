from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Alkene


class TerminalAlkeneFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.TerminalAlkene` instances.

    """

    _functional_group_smarts = '[*][C]([*])=[C]([H])[H]'
    _bonder_smarts = ['[$([C]=[C]([H])[H])]']
    _deleter_smarts = (
        ['[$([H][C]([H])=[C])]']*2 +
        ['[$([C](=[C])([H])[H])]']
    )

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Alkene(
                carbon1=atoms[1],
                atom1=atoms[0],
                atom2=atoms[2],
                carbon2=atoms[3],
                atom3=atoms[4],
                atom4=atoms[5],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
