import numpy as np


def get_centroid(position_matrix, atom_ids):
    return np.sum(position_matrix[atom_ids, :], axis=0) / len(atom_ids)
