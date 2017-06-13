from ..base import Topology


class BonderLocation:
    def __init__(self, coord):
        self.coord = coord


class Node:
    def __init__(self, coord, bonder_locations):
        self.coord = coord
        self.bonder_locations = [BonderLocation(x) for
                                 x in bonder_locations]

class Hexagonal(Topology):
    positions = [Node([], [])]

    def place_mols(self, macro_mol):
        ...
