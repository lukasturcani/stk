"""
Thiol Factory
=============

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Thiol


class ThiolFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Thiol` instances.

    Creates functional groups from substructures, which match the
    ``[*][S][H]`` functional group string.

    Examples
    --------
    You want to create a building block which has :class:`.Thiol`
    functional groups. You want the sulfur atom in those functional
    groups to be the *bonder* atom, and the hydrogen atom to be the
    *deleter* atom.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='SCCCCS',
            functional_groups=(stk.ThiolFactory(), ),
        )

    You want to create a building block which has :class:`.Thiol`
    functional groups. You want the non-hydrogen atom bonded to the
    sulfur to be the *bonder* atom and the SH group to be *deleter*
    atoms.

    .. code-block:: python

        import stk

        thiol_factory = stk.ThiolFactory(
            # The index of the atom bonded to sulfur is 0 in the
            # functional group string (see docstring).
            bonders=(0, ),
            # The index of the sulfur and hydrogen atoms in the
            # functional group string (see docstring) is 1 and 2,
            # respectively.
            deleters=(1, 2),
        )
        building_block = stk.BuildingBlock(
            smiles='SCCCCS',
            functional_groups=(thiol_factory, ),
        )


    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(1, ), deleters=(2, )):
        """
        Initialize a :class:`.ThiolFactory` instance.

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

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids('[*][S][H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Thiol(
                sulfur=atoms[1],
                hydrogen=atoms[2],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
