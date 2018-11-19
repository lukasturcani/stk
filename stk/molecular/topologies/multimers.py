"""

"""

from .base import Topology
from ...utilities import translation_component, quaternion
import numpy as np
import rdkit.Chem.AllChem as rdkit


class Dimer(Topology):
    def __init__(self, u, distance):
        """

        """

        self.u = u
        self.distance = distance

    def place_mols(self, macro_mol):
        bb = macro_mol.building_blocks[0]
        bb.set_position([0, 0, 0])
        macro_mol.mol = rdkit.Mol(bb.mol)

        # The axis of rotation shoudl be the longest dimension of the
        # building block.
        _, a1, a2 = bb.max_diameter()
        a1_coord = bb.atom_coords(int(a1))
        a2_coord = bb.atom_coords(int(a2))
        axis = a1_coord - a2_coord
        bb.rotate(np.random.random()*np.pi/2, axis)
        q = quaternion(self.u)
        v = translation_component(q)*self.distance
        new_pos = (bb.mol.GetConformer().GetPositions() +
                   np.matrix(v)).T
        bb.set_position_from_matrix(new_pos)
        macro_mol.mol = rdkit.CombineMols(macro_mol.mol, bb.mol)

    def bonded_fgs(self, macro_mol):
        return []
