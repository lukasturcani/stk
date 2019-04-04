"""
Defines cage topologies of building blocks with 3 functional groups.

"""

import numpy as np
from scipy.spatial.distance import euclidean

from .base import NoLinkerCageTopology,  Vertex


class OnePlusOne(NoLinkerCageTopology):
    """
    A sandwich cage topology from tri-functionalised building blocks.

    """

    x = 1
    positions_A = [Vertex(x, 0., 0.),
                   Vertex(-x, 0., 0.)]
    a, b = positions_A
    connections = [(a, b)]

    a.edge_plane_normal = lambda scale: scale*np.array([1, 0, 0])
    b.edge_plane_normal = lambda scale: scale*np.array([-1, 0, 0])

    a.edge_centroid = lambda scale: scale*np.array([0, 0, 0])
    b.edge_centroid = lambda scale: scale*np.array([0, 0, 0])

    n_windows = 3
    n_window_types = 1

    def bonded_fgs(self, macro_mol):

        for position in self.positions_A:
            other_position = next(x for x in self.positions_A if
                                  x is not position)

            position.fg_position_pairs = [
                (fg, other_position) for fg in position.fgs
            ]

            for fg1, vertex in position.fg_position_pairs:
                # Get all the distances between the fg and the fgs
                # on the vertex. Store this information on the vertex.

                for fg2 in vertex.fgs:
                    c1 = macro_mol.atom_centroid(fg1.bonder_ids)
                    c2 = macro_mol.atom_centroid(fg2.bonder_ids)
                    distance = euclidean(c1, c2)
                    position.distances.append((distance, fg1, fg2))

        paired = set()
        for position in self.positions_A:
            for _, fg1, fg2 in sorted(position.distances):
                if fg1 in paired or fg2 in paired:
                    continue

                # Add the bond.
                yield fg1, fg2
                paired.add(fg1)
                paired.add(fg2)


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
