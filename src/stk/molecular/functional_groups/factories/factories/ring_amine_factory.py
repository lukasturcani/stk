from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import RingAmine


class RingAmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.RingAmine` functional groups.

    """

    def __init__(self):
        super().__init__(
            functional_group_smarts=(
                '[N]([H])([H])[#6]~[#6]([H])~[#6R1]'
            ),
            bonders=(),
            deleters=(),
        )

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield RingAmine(
                nitrogen=atoms[0],
                hydrogen1=atoms[1],
                hydrogen2=atoms[2],
                carbon1=atoms[3],
                carbon2=atoms[4],
                hydrogen3=atoms[5],
                carbon3=atoms[6],
            )
