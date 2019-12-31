from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Bromo


class BromoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Bromo` instances.

    """

    _functional_group_smarts = '*[Br]'
    _bonder_smarts = ['[$(*[Br])]']
    _deleter_smarts = ['[$([Br]*)]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Bromo(
                bromine=atoms[1],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )
