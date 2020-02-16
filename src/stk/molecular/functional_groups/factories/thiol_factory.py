from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Thiol


class ThiolFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Thiol` instances.

    Creates functional groups from substructures, which match the
    ``[*][S][H]`` functional group string.

    """

    def __init__(self, bonders=(1, ), deleters=(2, )):
        """
        Initialize a :class:`.ThiolFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are bonder atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are deleter atoms.

        """

        super().__init__('[*][S][H]', bonders, deleters)

    def get_functional_groups(self, molecule):
        for atom_ids in self._get_atom_ids(molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Thiol(
                sulfur=atoms[1],
                hydrogen=atoms[2],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
