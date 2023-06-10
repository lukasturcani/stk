"""
Dibromo Factory
===============

"""


from stk._internal.functional_group_factories.functional_group_factory import (
    FunctionalGroupFactory,
)
from stk._internal.functional_group_factories.utilities import get_atom_ids
from stk._internal.functional_groups.dibromo import Dibromo


class DibromoFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Dibromo` instances.

    Creates functional groups from substructures, which match the
    ``[Br][#6]~[#6][Br]`` functional group string.

    Examples
    --------
    *Creating Functional Groups with the Factory*

    You want to create a building block which has :class:`.Dibromo`
    functional groups. You want the non-bromine atoms of the functional
    group to be the *bonder* atoms and the bromine atoms to be the
    *deleter* atoms.

    .. testcode:: creating-functional-groups-with-the-factory

        import stk

        building_block = stk.BuildingBlock(
            smiles='BrCC(Br)CCCC',
            functional_groups=(stk.DibromoFactory(), ),
        )

    .. testcode:: creating-functional-groups-with-the-factory
        :hide:

        assert all(
            isinstance(functional_group, stk.Dibromo)
            for functional_group
            in building_block.get_functional_groups()
        )
        assert building_block.get_num_functional_groups() == 1

    *Changing the Bonder and Deleter Atoms*

    You want to create a building block which has :class:`.Dibromo`
    functional groups, You want only one of non-bromine atoms to be
    a *bonder* atom and its neighboring bromine atom to be a
    *deleter* atom.

    .. testcode:: changing-the-bonder-and-deleter-atoms

        import stk

        dibromo_factory = stk.DibromoFactory(
            # The index of one of the non-bromine atoms in the
            # functional group string (see docstring) is 1.
            bonders=(1, ),
            # The neighboring bromine atom has an index of 0.
            deleters=(0, ),
        )
        building_block = stk.BuildingBlock(
            smiles='BrCC(Br)CCC',
            functional_groups=(dibromo_factory, ),
        )

    .. testcode:: changing-the-bonder-and-deleter-atoms
        :hide:

        fg, = building_block.get_functional_groups()
        assert fg.get_num_bonders() == 1
        assert sum(1 for _ in fg.get_deleters()) == 1

        assert all(
            isinstance(atom, stk.C)
            for functional_group
            in building_block.get_functional_groups()
            for atom
            in functional_group.get_bonders()
        )
        assert all(
            isinstance(atom, stk.Br)
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

    def __init__(self, bonders=(1, 2), deleters=(0, 3), placers=None):
        """
        Initialize a :class:`.DibromoFactory` instance.

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

        self._bonders = bonders
        self._deleters = deleters
        self._placers = bonders if placers is None else placers

    def get_functional_groups(self, molecule):
        for atom_ids in get_atom_ids("[Br][#6]~[#6][Br]", molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Dibromo(
                atom1=atoms[1],
                bromine1=atoms[0],
                atom2=atoms[2],
                bromine2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
