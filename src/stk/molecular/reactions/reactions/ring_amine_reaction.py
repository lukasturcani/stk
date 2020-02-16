import numpy as np

from .reaction import Reaction
from ... import atoms
from ...bond import Bond


class RingAmineReaction(Reaction):
    def __init__(self, position_matrix, ring_amine1, ring_amine2):
        self._position_matrix = np.array(position_matrix)
        self._ring_amine1 = ring_amine1.clone()
        self._ring_amine2 = ring_amine2.clone()

    def _get_new_atoms(self):
        n1_coord = self._get_position(self._ring_amine1.get_nitrogen())
        n2_coord = self._get_position(self._ring_amine2.get_nitrogen())
        c1_coord = self._get_position(self._ring_amine1.get_carbon1())
        c2_coord = self._get_position(self._ring_amine2.get_carbon1())

        n_joiner = atoms.C(-1)
        n_joiner_coord = (n1_coord + n2_coord) / 2
        yield n_joiner, n_joiner_coord

        nh1 = atoms.H(-2)
        nh1_coord = n_joiner_coord + [0, 0, 1]
        yield nh1, nh1_coord

        nh2 = atoms.H(-3)
        nh2_coord = n_joiner_coord + [0, 0, -1]
        yield nh2, nh2_coord

        nc_joiner1 = atoms.C(-4)
        nc_joiner1_coord = (c1_coord + n2_coord) / 2
        yield nc_joiner1, nc_joiner1_coord

        nc1h1 = atoms.H(-5)
        nc1h1_coord = nc_joiner1_coord + [0, 0, 1]
        yield nc1h1, nc1h1_coord

        nc1h2 = atoms.H(-6)
        nc1h2_coord = nc_joiner1_coord + [0, 0, -1]
        yield nc1h2, nc1h2_coord

        nc_joiner2 = atoms.C(-7)
        nc_joiner2_coord = (c2_coord + n1_coord) / 2
        yield nc_joiner2, nc_joiner2_coord

        nc2h1 = atoms.H(-8)
        nc2h1_coord = nc_joiner2_coord + [0, 0, 1]
        yield nc2h1, nc2h1_coord

        nc2h2 = atoms.H(-9)
        nc2h2_coord = nc_joiner2_coord + [0, 0, -1]
        yield nc2h2, nc2h2_coord

    def _get_new_bonds(self):
        n1 = self._ring_amine1.get_nitrogen()
        n2 = self._ring_amine2.get_nitrogen()
        c1 = self._ring_amine1.get_carbon2()
        c2 = self._ring_amine2.get_carbon2()
        n_joiner = atoms.C(-1)
        nh1 = atoms.H(-2)
        nh2 = atoms.H(-3)
        nc_joiner1 = atoms.C(-4)
        nc1h1 = atoms.H(-5)
        nc1h2 = atoms.H(-6)
        nc_joiner2 = atoms.C(-7)
        nc2h1 = atoms.C(-8)
        nc2h2 = atoms.C(-9)

        yield Bond(n1, n_joiner, 1),
        yield Bond(n_joiner, n2, 1, self._periodicity),
        yield Bond(n_joiner, nh1, 1),
        yield Bond(n_joiner, nh2, 1),
        yield Bond(c1, nc_joiner1, 1),
        yield Bond(nc_joiner1, n2, 1, self._periodicity),
        yield Bond(nc_joiner1, nc1h1, 1),
        yield Bond(nc_joiner1, nc1h2, 1),
        yield Bond(nc_joiner2, c2, 1, self._periodicity),
        yield Bond(n1, nc_joiner2, 1),
        yield Bond(nc_joiner2, nc2h1, 1),
        yield Bond(nc_joiner2, nc2h2, 1)

    def _get_deleted_atoms(self):
        yield self._ring_amine1.get_hydrogen3()
        yield self._ring_amine2.get_hydrogen3()
        yield self._ring_amine1.get_hydrogen1()
        yield self._ring_amine2.get_hydrogen1()
        yield self._ring_amine1.get_hydrogen2()
        yield self._ring_amine2.get_hydrogen2()
