"""
Boronic Acid Factory
====================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import BoronicAcid


class BoronicAcidFactory(FunctionalGroupFactory):
    """
    Creates :class:`.BoronicAcid` instances.

    Creates functional groups from substructures, which match the
    ``[*][B]([O][H])[O][H]`` functional group string.

    Examples
    --------
    You want to create a building block which has :class:`.BoronicAcid`
    functional groups. You want the boron atom in those functional
    groups to be the *bonder* atom and the OH groups to be *deleter*
    atoms.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='OB(O)CCCB(O)O',
            functional_groups=(stk.BoronicAcidFactory(), ),
        )

    You want to create a building block which has :class:`.BoronicAcid`
    functional groups. You want the oxygen atoms to be treated as
    *bonder* atoms, and the hydrogen atoms to be treated as *deleter*
    atoms.

    .. code-block:: python

        import stk

        boronic_acid_factory = stk.BoronicAcidFactory(
            # The indices of the oxygen atoms in the functional
            # group string (see docstring) are 2 and 4.
            bonders=(2, 4),
            # The indices of the hydrogen atoms in the
            # functional group string (see docstring) are 3 and 5.
            deleters=(3, 5),
        )
        building_block = stk.BuildingBlock(
            smiles='OB(O)CCCB(O)O',
            functional_groups=(boronic_acid_factory, ),
        )


    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(
        self,
        bonders=(1, ),
        deleters=(2, 3, 4, 5),
        placers=None,
    ):
        """
        Initialize a :class:`.BoronicAcidFactory` instance.

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
        ids = _get_atom_ids('[*][B]([O][H])[O][H]', molecule)
        for atom_ids in ids:
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield BoronicAcid(
                boron=atoms[1],
                oxygen1=atoms[2],
                hydrogen1=atoms[3],
                oxygen2=atoms[4],
                hydrogen2=atoms[5],
                atom=atoms[0],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
