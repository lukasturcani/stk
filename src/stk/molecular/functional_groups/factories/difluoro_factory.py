"""
Difluoro Factory
================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Difluoro


class DifluoroFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Difluoro` instances.

    Creates functional groups from substructures, which match the
    ``[F][#6]~[#6][F]`` functional group string.

    Examples
    --------
    You want to create a building block which has :class:`.Difluoro`
    functional groups. You want the non-fluorine atoms in those
    functional groups to be the *bonder* atoms, and the fluorine atoms
    to be the *deleter* atoms.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='FCC(F)CCC',
            functional_groups=(stk.Difluoro(), ),
        )

    You want to create a building block which has :class:`.Difluoro`
    functional groups, You want only one of non-fluorine atoms to be
    a *bonder* atom and its neighboring fluorine atom to be a
    *deleter* atom.

    .. code-block:: python

        import stk

        difluoro_factory = stk.DifluoroFactory(
            # The index of one of the non-fluorine atoms in the
            # functional group string (see docstring) is 1.
            bonders=(1, ),
            # The neighboring fluorine atom has an index of 0.
            deleters=(0, ),
        )
        building_block = stk.BuildingBlock(
            smiles='FCC(F)CCC',
            functional_groups=(difluoro_factory, ),
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(1, 2), deleters=(0, 3), placers=None):
        """
        Initialize a :class:`.DifluoroFactory` instance.

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
        for atom_ids in _get_atom_ids('[F][#6]~[#6][F]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Difluoro(
                atom1=atoms[1],
                fluorine1=atoms[0],
                atom2=atoms[2],
                fluorine2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
