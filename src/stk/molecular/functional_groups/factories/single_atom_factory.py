"""
Single Atom Factory
===================
"""

from .functional_group_factory import FunctionalGroupFactory
from ..functional_groups import SingleAtom


class SingleAtomFactory(FunctionalGroupFactory):
    """
    Creates :class:`.SingleAtom` instances.

    This factory does not use predefined substructures to find the
    functional groups. Instead it supplies a desired number of
    functional groups for a single atom.

    Examples
    --------
    You want to create a single atom building block which has
    :class:`.SingleAtom` functional groups.

    .. code-block:: python
        import stk

        # Build single atoms to place on metal centre topology.
        # Metal atom with 6 functional groups.
        atom = rdkit.MolFromSmiles('[Fe+2]')
        atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
        metal_atom = stk.BuildingBlock.init_from_rdkit_mol(
            atom,
            functional_groups=[stk.SingleAtomFactory(
                num_functional_groups=6
            )]
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, num_functional_groups=1):
        """
        Initialize a :class:`.SingleAtomFactory` instance.

        Parameters
        ----------
        num_functional_groups : :class:`int`
            The target number of functional groups to apply to the
            single atom.

        """

        self._num_functional_groups = num_functional_groups

    def get_functional_groups(self, molecule):
        for i in range(self._num_functional_groups):
            atoms = tuple(molecule.get_atoms((0)))
            yield SingleAtom(atom=atoms[0])
