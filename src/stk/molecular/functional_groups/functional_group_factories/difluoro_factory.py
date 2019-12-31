from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Difluoro


class DifluoroFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Difluoro` instances.

    """

    _functional_group_smarts = '[F][#6]~[#6][F]'
    _bonder_smarts = ['[$([#6]([F])~[#6][F])]']*2
    _deleter_smarts = ['[$([F][#6]~[#6][F])]']*2

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Difluoro(
                atom1=atoms[1],
                fluorine1=atoms[0],
                atom2=atoms[2],
                fluorine2=atoms[3],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
