"""
Ring Amine Factory
==================

"""


from stk._internal.functional_group_factories.functional_group_factory import (
    FunctionalGroupFactory,
)
from stk._internal.functional_group_factories.utilities import get_atom_ids
from stk._internal.functional_groups.ring_amine import RingAmine


class RingAmineFactory(FunctionalGroupFactory):
    """
    Creates :class:`.RingAmine` functional groups.

    Creates functional groups from substructures, which match the
    ``[N]([H])([H])[#6]~[#6]([H])~[#6R1]`` functional group string.

    """

    # Keep this __init__() method, which does nothing, here for the
    # docstring.
    def __init__(self):
        """
        Initialize a :class:`RingAmineFactory` instance.

        """

        return super().__init__()

    def get_functional_groups(self, molecule):
        ids = get_atom_ids(
            query="[N]([H])([H])[#6]~[#6]([H])~[#6R1]",
            molecule=molecule,
        )
        for atom_ids in ids:
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
