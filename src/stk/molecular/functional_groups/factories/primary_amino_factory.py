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
    You want to create a building block which has
    :class:`.PrimaryAmino` functional groups. You want the nitrogen
    atom to be the *bonder* atom, and the hydrogen atoms to be the
    *deleter* atoms.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='NCCCCN',
            functional_groups=(stk.PrimaryAminoFactory(), ),
        )

    You want to create a building block which has
    :class:`.PrimaryAmino` functional groups. You want the non-hydrogen
    atom bonded to nitrogen to be the *bonder* atom and the
    nitrogen and hydrogen atoms to be *deleter* atoms.

    .. code-block:: python

        import stk

        primary_amino_factory = stk.PrimaryAminoFactory(
            # The index of the atom attached to the nitrogen is 0 in
            # the functional group string (see docstring).
            bonders=(0, ),
            # The indices of the nitrogen and hydrogen atoms in the
            # functional group string (see docstring) are 1, 2 and 3.
            deleters=(1, 2, 3),
        )
        building_block = stk.BuildingBlock(
            smiles='NCCCCN',
            functional_groups=(primary_amino_factory, ),
        )


    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3), placers=None):
        """
        Initialize a :class:`.PrimaryAminoFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are *bonder* atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are *deleter* atoms.

        placers : :class:`tuple` of :class:`int`, optional
            The indices of atoms in the functional group string, which
            are *placer* atoms. If ``None``, `bonders` will be used.

        """

        """
        Initialize an :class:`.AmineFactory`.

        """

        self._bonders = bonders
        self._deleters = deleters
        self._placers = bonders if placers is None else placers

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
                placers=tuple(atoms[i] for i in self._placers),
            )
