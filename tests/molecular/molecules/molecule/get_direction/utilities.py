import numpy as np


def get_direction(position_matrix, atom_ids):
    atom_positions = position_matrix[atom_ids, :]
    centered_positions = atom_positions - atom_positions.mean(axis=0)
    return np.linalg.svd(centered_positions)[-1][0]
