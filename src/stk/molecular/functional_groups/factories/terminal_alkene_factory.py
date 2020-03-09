"""
Terminal Alkene Factory
=======================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Alkene


class TerminalAlkeneFactory(FunctionalGroupFactory):
    """
    Creates :class:`.TerminalAlkene` instances.

    Creates functional groups from substructures, which match the
    ``[*][C]([*])=[C]([H])[H]`` functional group string.

    Examples
    --------
    You want to create a building block which has
    :class:`.TerminalAlkene` functional groups. You want the
    non-terminal carbon atom in those functional
    groups to be the *bonder* atom, and the terminal CH:subscript:`2`
    group to be the *deleter* atoms.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='C=CCCCCC=C',
            functional_groups=(stk.TerminalAlkeneFactory(), ),
        )

    You want to create a building block which has
    :class:`.TerminalAlkene` functional groups. You want the carbon
    atoms to be the *bonder* atoms and you don't want any *deleter*
    atoms.

    .. code-block:: python

        import stk

        terminal_alkene_factory = stk.TerminalAlkeneFactory(
            # The index of the carbon atoms in the functional
            # group string (see docstring) are 1 and 3.
            bonders=(1, 3),
            deleters=(),
        )
        building_block = stk.BuildingBlock(
            smiles='C=CCCCCC=C',
            functional_groups=(terminal_alkene_factory, ),
        )


    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, bonders=(1, ), deleters=(3, 4, 5)):
        """
        Initialize a :class:`.TerminalAlkeneFactory` instance.

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
        ids = _get_atom_ids('[*][C]([*])=[C]([H])[H]', molecule)
        for atom_ids in ids:
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Alkene(
                carbon1=atoms[1],
                atom1=atoms[0],
                atom2=atoms[2],
                carbon2=atoms[3],
                atom3=atoms[4],
                atom4=atoms[5],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )
