from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Diol


class DiolFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Diol` instances.

    """

    _functional_group_smarts = '[H][O][#6]~[#6][O][H]'
    _bonder_smarts = ['[$([O]([H])[#6]~[#6][O][H])]']*2
    _deleter_smarts = ['[$([H][O][#6]~[#6][O][H])]']*2

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Diol(
                hydrogen1=atoms[0],
                oxygen1=atoms[1],
                atom1=atoms[2],
                atom2=atoms[3],
                oxygen2=atoms[4],
                hydrogen2=atoms[5],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
