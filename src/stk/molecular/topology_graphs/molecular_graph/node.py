import numpy as np


class Node:
    def __init__(self, building_block, position=(0., 0., 0.)):
        self._building_block = building_block
        self._position = np.array(position)

    def get_building_block(self):
        return self._building_block

    def get_position(self):
        return np.array(self._position)
