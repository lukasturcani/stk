"""
Defines cage topologies of building blocks with 3 functional groups.

"""

import numpy as np
import rdkit.Chem.AllChem as rdkit

from .base import NoLinkerCageTopology,  Vertex


class OnePlusOne(NoLinkerCageTopology):
    """
    A sandwich cage topology from tri-functionalised building blocks.

    """

    x = 1
    positions_A = [Vertex(x, 0, 0),
                   Vertex(-x, 0, 0)]
    a, b = positions_A
    connections = [(a, b)]

    a.edge_plane_normal = lambda scale: scale*np.array([1, 0, 0])
    b.edge_plane_normal = lambda scale: scale*np.array([-1, 0, 0])

    a.edge_centroid = lambda scale: scale*np.array([0, 0, 0])
    b.edge_centroid = lambda scale: scale*np.array([0, 0, 0])

    n_windows = 3
    n_window_types = 1

    def join_mols(self, macro_mol):

        editable_mol = rdkit.EditableMol(macro_mol.mol)

        for position in self.positions_A:
            other_position = next(x for x in self.positions_A if
                                  x is not position)

            position.atom_position_pairs = [(atom, other_position) for
                                            atom in position.bonder_ids]

            for atom_id, vertex in position.atom_position_pairs:
                # Get all the distances between the atom and the bonder
                # atoms on the vertex. Store this information on the
                # vertex.
                for atom2_id in vertex.bonder_ids:
                    distance = macro_mol.atom_distance(atom_id,
                                                       atom2_id)
                    position.distances.append((distance,
                                              atom_id, atom2_id))

        paired = set()
        for position in self.positions_A:
            for _, atom1_id, atom2_id in sorted(position.distances):
                if atom1_id in paired or atom2_id in paired:
                    continue

                bond_type = self.determine_bond_type(macro_mol,
                                                     atom1_id,
                                                     atom2_id)
                # Add the bond.
                editable_mol.AddBond(atom1_id, atom2_id, bond_type)
                macro_mol.bonds_made += 1
                paired.add(atom1_id)
                paired.add(atom2_id)

        macro_mol.mol = editable_mol.GetMol()


class TwoPlusTwo(NoLinkerCageTopology):
    """
    Tetrahedral cage topology from tri-functionalised building blocks.

    """

    x = 1
    positions_A = [Vertex(x, 0, -x/np.sqrt(2)),
                   Vertex(-x, 0, -x/np.sqrt(2)),
                   Vertex(0, x, x/np.sqrt(2)),
                   Vertex(0, -x, x/np.sqrt(2))]

    a, b, c, d = positions_A

    for x in positions_A:
        old_normal = x.edge_plane_normal
        x.edge_plane_normal = lambda scale, a=old_normal: -1*a(scale)

    connections = [(a, b), (a, c), (a, d),
                   (b, c), (b, d),
                   (c, d)]

    n_windows = 4
    n_window_types = 1


class FourPlusFour(NoLinkerCageTopology):
    """
    A square cage topology from tri-functionalised building blocks.

    """

    x = 1
    positions_A = [Vertex(-x, x, -x),
                   Vertex(-x, -x, -x),
                   Vertex(x, x, -x),
                   Vertex(x, -x, -x),

                   Vertex(-x, x, x),
                   Vertex(-x, -x, x),
                   Vertex(x, x, x),
                   Vertex(x, -x, x)]

    a, b, c, d, e, f, g, h = positions_A

    connections = [(a, b), (a, c), (a, e), (b, d), (b, f), (c, g),
                   (c, d), (d, h), (e, g), (e, f), (f, h), (g, h)]

    n_windows = 6
    n_window_types = 1
