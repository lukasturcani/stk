from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Fluoro


class FluoroFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Fluoro` instances.

    """

    _functional_group_smarts = '*[F]'
    _bonder_smarts = ['[$(*[F])]']
    _deleter_smarts = ['[$([F]*)]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Fluoro(
                fluorine=atoms[1],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )
