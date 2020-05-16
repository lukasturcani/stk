"""
Diol Factory
============

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Diol


class DiolFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Diol` instances.

    Creates functional groups from substructures, which match the
    ``[H][O][#6]~[#6][O][H]`` functional group string.

    Examples
    --------
    You want to create a building block which has :class:`.Diol`
    functional groups. You want the carbon atoms in those functional
    groups to be the *bonder* atoms, and the OH groups to be a leaving
    groups.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='CCCC(O)C(O)CCCC',
            functional_groups=(stk.DiolFactory(), ),
        )

    You want to create a building block which has :class:`.Diol`
    functional groups. You want the oxygen atoms to be the *bonder*
    atoms and the hydrogen atoms to be the *deleter* atoms.

    .. code-block:: python

        import stk

        diol_factory = stk.DiolFactory(
            # The indices of the oxygen atoms in the functional
            # group string (see docstring) are 1 and 4.
            bonders=(1, 4),
            # The indices of the hydrogen atoms in the functional
            # group string (see docstring) are 0 and 5.
            deleters=(0, 5),
        )
        building_block = stk.BuildingBlock(
            smiles='CCCC(O)C(O)CCCC',
            functional_groups=(diol_factory, ),
        )


    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(
        self,
        bonders=(2, 3),
        deleters=(0, 1, 4, 5),
        placers=None,
    ):
        """
        Initialize a :class:`.DiolFactory` instance.

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
        ids = _get_atom_ids('[H][O][#6]~[#6][O][H]', molecule)
        for atom_ids in ids:
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Diol(
                hydrogen1=atoms[0],
                oxygen1=atoms[1],
                atom1=atoms[2],
                atom2=atoms[3],
                oxygen2=atoms[4],
                hydrogen2=atoms[5],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
