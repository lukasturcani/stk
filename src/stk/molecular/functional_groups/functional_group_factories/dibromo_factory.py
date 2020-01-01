from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Dibromo


class DibromoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Dibromo` instances.

    """

    _functional_group_smarts = '[Br][#6]~[#6][Br]'
    _bonder_smarts = ['[$([#6]([Br])~[#6][Br])]']*2
    _deleter_smarts = ['[$([Br][#6]~[#6][Br])]']*2

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Dibromo(
                atom1=atoms[1],
                bromine1=atoms[0],
                atom2=atoms[2],
                bromine2=atoms[3],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
