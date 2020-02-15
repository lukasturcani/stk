from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Alkene


class TerminalAlkeneFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.TerminalAlkene` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(3, 4, 5)):
        super().__init__('[*][C]([*])=[C]([H])[H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Alkene(
                carbon1=atoms[1],
                atom1=atoms[0],
                atom2=atoms[2],
                carbon2=atoms[3],
                atom3=atoms[4],
                atom4=atoms[5],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
