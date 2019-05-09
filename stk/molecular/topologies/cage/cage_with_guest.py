import rdkit.Chem.AllChem as rdkit

from ..base import Topology
from ....utilities import (rotation_matrix_arbitrary_axis,
                           add_fragment_props)


class CageWithGuest(Topology):
    """
    A topology representing a cage with a guest.

    When using this topology the cage must be the first building block
    in :attr:`.MacroMolecule.building_blocks` and the guest must be
    the second.

    Attributes
    ----------
    rot_mat : :class:`numpy.ndarray`
        The rotation matrix applied to the guest molecule.

    displacement : :class:`list` of :class:`float`
        The translational offset of the guest from the center of the
        cage cavity.

    """

    def __init__(self, axis=None, angle=0, displacement=None):
        """
        Initializes an instance of :class:`CageWithGuest`.

        Parameters
        ----------
        axis : :class:`list` of :class:`int`, optional
            The axis about which the guest is rotated.

        angle : :class:`float`
            The angle by which the guest is rotated.

        displacement : :class:`list` of :class:`float`
            The translational offset of the guest from the center of
            the cage cavity.

        """

        self.rot_mat = None
        if axis is not None:
            self.rot_mat = rotation_matrix_arbitrary_axis(angle, axis)

        self.displacement = None
        if displacement is not None:
            # Nested list so that it can be added to a [3, n] array.
            self.displacement = [[i] for i in displacement]

        super().__init__(del_atoms=False)

    def place_mols(self, macro_mol):
        origin = [0, 0, 0]

        for i, bb in enumerate(macro_mol.building_blocks):
            pos = bb.position_matrix()
            bb.set_position(origin)

            # Make sure that the building blocks are always positioned
            # consistently, regardless of initial position.
            _, a1, a2 = bb.max_diameter()
            mol_axis = bb.atom_coords(a1) - bb.atom_coords(a2)
            bb.set_orientation(mol_axis, [0, 1, 0])

            # Apply the rotation and displacement to the guest.
            if i == 1:
                new_pos = bb.position_matrix()
                if self.rot_mat is not None:
                    new_pos = self.rot_mat @ new_pos
                if self.displacement is not None:
                    new_pos += self.displacement

                bb.set_position_from_matrix(new_pos)

            # Make a copy so that atoms can be tagged without modifying
            # the original.
            bb_mol = rdkit.Mol(bb.mol)
            add_fragment_props(bb_mol, i, i)
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol, bb_mol)
            bb.set_position_from_matrix(pos)

    def bonded_fgs(self, macro_mol):
        return iter(())
