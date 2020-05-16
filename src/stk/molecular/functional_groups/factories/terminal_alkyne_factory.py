"""
Terminal Alkyne Factory
=======================

"""
from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Alkyne


class TerminalAlkyneFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Alkyne` instances.

    Creates functional groups from substructures, which match the
    ``[*][C]#[C][H]`` functional group string.

    Examples
    --------
    You want to create a building block which has
    :class:`.Alkyne` functional groups, but only if they are terminal.
    You want the non-terminal carbon atom in those functional
    groups to be the *bonder* atom, and the terminal CH
    group to be the *deleter* atoms.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='C#CCCCCC#C',
            functional_groups=(stk.TerminalAlkyneFactory(), ),
        )

    You want to create a building block which has
    :class:`.Alkyne` functional groups. You want the carbon
    atoms to be the *bonder* atoms and you don't want any *deleter*
    atoms.

    .. code-block:: python

        import stk

        terminal_alkyne_factory = stk.TerminalAlkyneFactory(
            # The indices of the carbon atoms in the functional
            # group string (see docstring) are 1 and 2.
            bonders=(1, 2),
            deleters=(),
        )
        building_block = stk.BuildingBlock(
            smiles='C#CCCCCC#C',
            functional_groups=(terminal_alkyne_factory, ),
        )


    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(1, ), deleters=(2, 3), placers=None):
        """
        Initialize a :class:`.TerminalAlkyneFactory` instance.

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
        for atom_ids in _get_atom_ids('[*][C]#[C][H]', molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Alkyne(
                atom1=atoms[0],
                carbon1=atoms[1],
                carbon2=atoms[2],
                atom2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
