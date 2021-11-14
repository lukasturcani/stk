"""
Secondary Amino Factory
=======================

"""

from ..functional_groups import SecondaryAmino
from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids


class SecondaryAminoFactory(FunctionalGroupFactory):
    """
    Creates :class:`.SecondaryAmino` instances.

    Creates functional groups from substructures, which match the
    ``[H][N]([#6])[#6]`` functional group string.

    Examples
    --------
    *Creating Functional Groups with the Factory*

    You want to create a building block which has
    :class:`.SecondaryAmino` functional groups. You want the nitrogen
    atom in those functional groups to be the *bonder* atom, and the
    hydrogen atom to be the *deleter* atom.

    .. testcode:: creating-functional-groups-with-the-factory

        import stk

        building_block = stk.BuildingBlock(
            smiles='N(CCCC)CCCC',
            functional_groups=(stk.SecondaryAminoFactory(), ),
        )

    .. testcode:: creating-functional-groups-with-the-factory
        :hide:

        assert all(
            isinstance(functional_group, stk.SecondaryAmino)
            for functional_group
            in building_block.get_functional_groups()
        )
        assert building_block.get_num_functional_groups() == 1

    *Changing the Bonder and Deleter Atoms*

    You want to create a building block which has
    :class:`.SecondaryAmino` functional groups. You want the nitrogen
    atom to be the *bonder* atom and one of the carbon atoms to be the
    *deleter* atom.

    .. testcode:: changing-the-bonder-and-deleter-atoms

        import stk

        secondary_amino_factory = stk.SecondaryAminoFactory(
            # The index of the nitrogen atom in the functional
            # group string (see docstring) is 1.
            bonders=(1, ),
            # The index of one of the carbon atoms in the functional
            # group string (see docstring) is 2.
            deleters=(2, ),
        )
        building_block = stk.BuildingBlock(
            smiles='N(CCCC)CCCC',
            functional_groups=(secondary_amino_factory, ),
        )

    .. testcode:: changing-the-bonder-and-deleter-atoms
        :hide:

        fg1, = building_block.get_functional_groups()
        assert fg1.get_num_bonders() == 1
        assert sum(1 for _ in fg1.get_deleters()) == 1

        assert all(
            isinstance(atom, stk.N)
            for functional_group
            in building_block.get_functional_groups()
            for atom
            in functional_group.get_bonders()
        )
        assert all(
            isinstance(atom, stk.C)
            for functional_group
            in building_block.get_functional_groups()
            for atom
            in functional_group.get_deleters()
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(1, ), deleters=(0, ), placers=None):
        """
        Initialize a :class:`.SecondaryAminoFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are bonder atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are deleter atoms.

        placers : :class:`tuple` of :class:`int`, optional
            The indices of atoms in the functional group string, which
            are *placer* atoms. If ``None``, `bonders` will be used.

        """

        self._bonders = bonders
        self._deleters = deleters
        self._placers = bonders if placers is None else placers

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids('[H][N]([#6])[#6]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield SecondaryAmino(
                nitrogen=atoms[1],
                hydrogen=atoms[0],
                atom1=atoms[2],
                atom2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
