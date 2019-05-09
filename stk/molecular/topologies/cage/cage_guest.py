from ..base import Topology
from ....utilities import rotation_matrix_arbitrary_axis


class CageGuest(Topology):
    """

    """

    def __init__(self, axis=None, angle=0, displacement=None):
        """

        """

        self.rot_mat = None
        if axis is not None:
            self.rot_mat = rotation_matrix_arbitrary_axis(angle, axis)
        
        self.displacement = None
        if displacement is not None:
            # Nested list so that it can be added to a [3, n] array.
            self.displacement = [[i] for i in displacement]

    def place_mols(self, macro_mol):
        origin = [0, 0, 0]

        for i, bb in enumerate(macro_mol.building_blocks):
            pos = bb.position_matrix()
            bb.set_position(origin)
            
            # Make sure that the building blocks are always positioned
            # consistently, regardless of initial position.
            _, a1, a2 = bb.max_diameter()
            mol_axis = bb.atom_coords(a1) - bb.atom_coords(a2)
            bb.set_orientation(mol_axis, [1, 0, 0])

            # Apply the rotation and displacement to the guest.
            if i == 1:
                new_pos = bb.position_matrix()
                if self.rot_mat is not None:
                    new_pos = self.rot_mat @ new_pos
                if self.displacement is not None:
                    new_pos += self.displacement

                bb.set_position_from_matrix(new_pos)

            macro_mol.mol = rdkit.CombineMols(macro_mol.mol, bb.mol)
            bb.set_position_from_matrix(pos)

    def bonded_fgs(self, macro_mol):
        return iter(())      