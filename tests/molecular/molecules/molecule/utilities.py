import numpy as np


def get_centroid(position_matrix, atom_ids):
    return np.sum(position_matrix[atom_ids, :], axis=0) / len(atom_ids)


def get_direction(position_matrix, atom_ids):
    atom_positions = position_matrix[atom_ids, :]
    centered_positions = atom_positions - atom_positions.mean(axis=0)
    return np.linalg.svd(centered_positions)[-1][0]
