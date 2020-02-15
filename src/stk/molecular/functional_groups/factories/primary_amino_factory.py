from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import PrimaryAmino


class PrimaryAminoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.PrimaryAmino` instances.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3)):
        """
        Initialize an :class:`.AmineFactory`.

        """

        super().__init__('[*][N]([H])[H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield PrimaryAmino(
                nitrogen=atoms[1],
                hydrogen1=atoms[2],
                hydrogen2=atoms[3],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
