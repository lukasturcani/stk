"""
SMARTS Functional Group Factory
===============================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import GenericFunctionalGroup


class SmartsFunctionalGroupFactory(FunctionalGroupFactory):
    """
    Creates :class:`.GenericFunctionalGroup` instances.

    Examples
    --------
    *Using SMARTS to Define Functional Groups*

    You want to create a building block which has
    :class:`.GenericFunctionalGroup` functional groups based on the
    SMARTS string: ``[Br][C]``.
    You want the ``C`` atom to be the *bonder* atom, and the
    ``Br`` atom to be the *deleter* atom.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='BrCCCBr',
            functional_groups=(
                stk.SmartsFunctionalGroupFactory(
                    smarts='[Br][C]',
                    bonders=(1, ),
                    deleters=(0, ),
                ),
            ),
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, smarts, bonders, deleters, placers=None):
        """
        Initialize a :class:`.SmartsFunctionalGroupFactory` instance.

        Parameters
        ----------
        smarts : :class:`str`
            The SMARTS defining the functional group.

        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in `smarts`, which are *bonder* atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in `smarts`, which are *deleter*
            atoms.

        placers : :class:`tuple` of :class:`int`, optional
            The indices of atoms in `smarts`, which are *placer* atoms.
            If ``None``, the *bonder* atoms will be used.

        """

        self._smarts = smarts
        self._bonders = bonders
        self._deleters = deleters
        self._placers = bonders if placers is None else placers

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids(self._smarts, molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield GenericFunctionalGroup(
                atoms=atoms,
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
