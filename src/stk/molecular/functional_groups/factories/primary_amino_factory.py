"""
Primary Amino Factory
=====================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import PrimaryAmino


class PrimaryAminoFactory(FunctionalGroupFactory):
    """
    Creates :class:`.PrimaryAmino` instances.

    Creates functional groups from substructures, which match the
    ``[*][N]([H])[H]`` functional group string.

    Examples
    --------
    You want to create a building block which has :class:`.Amide`
    functional groups. You want the carbon atom in those functional
    groups to be the bonder atom, and the amino group to be a leaving
    group.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='NC(=O)CC(=O)N',
            functional_groups=(stk.AmideFactory(), ),
        )

    You want to create a building block which has :class:`.Amide`
    functional groups. You want the carbon atom to be the bonder
    atom and the oxygen atom to be the deleter atom.

    .. code-block:: python

        import stk

        amide_factory = stk.AmideFactory(
            # The index of the carbon atom in the functional
            # group string (see docstring) is 1.
            bonders=(1, ),
            # The index of the oxygen atom in the functional
            # group string (see docstring) is 2.
            deleters=(2, ),
        )
        building_block = stk.BuildingBlock(
            smiles='NC(=O)CC(=O)N',
            functional_groups=(amide_factory, ),
        )


    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3)):
        """
        Initialize a :class:`.PrimaryAminoFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are bonder atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are deleter atoms.

        """

        """
        Initialize an :class:`.AmineFactory`.

        """

        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids('[*][N]([H])[H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield PrimaryAmino(
                nitrogen=atoms[1],
                hydrogen1=atoms[2],
                hydrogen2=atoms[3],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
