from ..base import Topology
from ....utilities import 


class CageGuest(Topoology):
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
            self.displacement = [[i] for i in displacement]

    def place_mols(self, macro_mol):
        origin = [0, 0, 0]

        for i, bb in enumerate(macro_mol.building_blocks):
            pos = bb.position_matrix()
            bb.set_position(origin)

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