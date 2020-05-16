"""
Fluoro Factory
==============

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Fluoro


class FluoroFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Fluoro` instances.

    Creates functional groups from substructures, which match the
    ``[*][F]`` functional group string.

    Examples
    --------
    You want to create a building block which has :class:`.Fluoro`
    functional groups. You want the non-fluorine atom in those
    functional groups to be the *bonder* atom, and the fluorine atom
    to be the *deleter* atom.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='FCCCF',
            functional_groups=(stk.FluoroFactory(), ),
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(0, ), deleters=(1, ), placers=None):
        """
        Initialize a :class:`.FluoroFactory` instance.

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
        for atom_ids in _get_atom_ids('[*][F]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Fluoro(
                fluorine=atoms[1],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
