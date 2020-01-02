from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ...functional_groups import Alkyne


class TerminalAlkyneFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Alkyne` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3)):
        super().__init__('[*][C]#[C][H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Alkyne(
                atom1=atoms[0],
                carbon1=atoms[1],
                carbon2=atoms[2],
                atom2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
