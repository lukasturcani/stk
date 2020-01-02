from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ...functional_groups import Bromo


class BromoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Bromo` instances.

    """

    def __init__(self, bonders=(0, ), deleters=(1, )):
        super().__init__('[*][Br]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Bromo(
                bromine=atoms[1],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
