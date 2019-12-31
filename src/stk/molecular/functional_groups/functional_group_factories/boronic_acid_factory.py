from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import BoronicAcid


class BoronicAcidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.BoronicAcid` instances.

    """

    _functional_group_smarts = '[*][B]([O][H])[O][H]'
    _bonder_smarts = ['[$([B]([O][H])[O][H])]']
    _deleter_smarts = (
        ['[$([O]([H])[B][O][H])]']*2 +
        ['[$([H][O][B][O][H])]']*2
    )

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield BoronicAcid(
                boron=atoms[1],
                oxygen1=atoms[2],
                hydrogen1=atoms[3],
                oxygen2=atoms[4],
                hydrogen2=atoms[5],
                atom=atoms[0],
                bonders=tuple(self._get_atoms(atoms, ids.bonder_ids)),
                deleters=tuple(
                    self._get_atoms(atoms, ids.deleter_ids)
                ),
            )
