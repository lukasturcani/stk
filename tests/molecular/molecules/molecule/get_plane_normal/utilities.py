import numpy as np


def get_plane_normal(position_matrix, atom_ids):
    atomic_positions = position_matrix[atom_ids, :]
    centroid = get_centroid(position_matrix, atom_ids)
    return np.linalg.svd(atomic_positions - centroid)[-1][2, :]


def get_centroid(position_matrix, atom_ids):
    return np.sum(position_matrix[atom_ids, :], axis=0) / len(atom_ids)
