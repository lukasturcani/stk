"""
Thioacid Factory
================

"""

from ..functional_groups import Thioacid
from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids


class ThioacidFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Thioacid` instances.

    Creates functional groups from substructures, which match the
    ``[*][C](=[O])[S][H]`` functional group string.

    Examples
    --------
    *Creating Functional Groups with the Factory*

    You want to create a building block which has :class:`.Thioacid`
    functional groups. You want the carbon atom in those functional
    groups to be the *bonder* atom, and SH group to be the *deleter*
    atoms.

    .. testcode:: creating-functional-groups-with-the-factory

        import stk

        building_block = stk.BuildingBlock(
            smiles='SC(=O)CC(=O)S',
            functional_groups=(stk.ThioacidFactory(), ),
        )

    .. testcode:: creating-functional-groups-with-the-factory
        :hide:

        assert all(
            isinstance(functional_group, stk.Thioacid)
            for functional_group
            in building_block.get_functional_groups()
        )
        assert building_block.get_num_functional_groups() == 2

    *Changing the Bonder and Deleter Atoms*

    You want to create a building block which has :class:`.Thioacid`
    functional groups. You want the carbon atom to be the *bonder*
    atom and the oxygen atom to be the *deleter* atom.

    .. testcode:: changing-the-bonder-and-deleter-atoms

        import stk

        thioacid_factory = stk.ThioacidFactory(
            # The index of the carbon atom in the functional
            # group string (see docstring) is 1.
            bonders=(1, ),
            # The index of the oxygen atom in the functional
            # group string (see docstring) is 2.
            deleters=(2, ),
        )
        building_block = stk.BuildingBlock(
            smiles='SC(=O)CC(=O)S',
            functional_groups=(thioacid_factory, ),
        )

    .. testcode:: changing-the-bonder-and-deleter-atoms
        :hide:

        fg1, fg2 = building_block.get_functional_groups()
        assert fg1.get_num_bonders() == 1
        assert sum(1 for _ in fg1.get_deleters()) == 1
        assert fg2.get_num_bonders() == 1
        assert sum(1 for _ in fg2.get_deleters()) == 1

        assert all(
            isinstance(atom, stk.C)
            for functional_group
            in building_block.get_functional_groups()
            for atom
            in functional_group.get_bonders()
        )
        assert all(
            isinstance(atom, stk.O)
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

    def __init__(self, bonders=(1, ), deleters=(3, 4), placers=None):
        """
        Initialize a :class:`.ThioacidFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are *bonder* atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are *deleter* atoms.

        """

        self._bonders = bonders
        self._deleters = deleters
        self._placers = bonders if placers is None else placers

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids('[*][C](=[O])[S][H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Thioacid(
                carbon=atoms[1],
                oxygen=atoms[2],
                sulfur=atoms[3],
                hydrogen=atoms[4],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
