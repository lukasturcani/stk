"""
Aldehyde Factory
================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Aldehyde


class AldehydeFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Aldehyde` instances.

    Creates functional groups from substructures, which match the
    ``[*][C](=[O])[H]`` functional group string.

    Examples
    --------
    You want to create a building block which has :class:`.Aldehyde`
    functional groups. You want the carbon atom in those functional
    groups to be the bonder atom, and the oxygen atom to be the
    deleter atom.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='O=CCC=O',
            functional_groups=(stk.AldehydeFactory(), ),
        )

    You want to create a building block which has :class:`.Aldehyde`
    functional groups. You want the carbon atom to be the bonder atom
    and the hydrogen atom to be the deleter atom.

    .. code-block:: python

        import stk

        aldehyde_factory = stk.AldehydeFactory(
            # The index of the carbon atom in the functional
            # group string (see docstring) is 1.
            bonders=(1, ),
            # The index of the hydrogen atom in the functional
            # group string (see docstring) is 3.
            deleters=(3, ),
        )
        building_block = stk.BuildingBlock(
            smiles='O=CCC=O',
            functional_groups=(aldehyde_factory, ),
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(1, ), deleters=(2, ), placers=None):
        """
        Initialize a :class:`.AldehydeFactory` instance.

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
        for atom_ids in _get_atom_ids('[*][C](=[O])[H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Aldehyde(
                carbon=atoms[1],
                oxygen=atoms[2],
                hydrogen=atoms[3],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
