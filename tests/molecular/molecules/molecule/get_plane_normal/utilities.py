import numpy as np
from ..utilities import get_centroid


def get_plane_normal(position_matrix, atom_ids):
    atomic_positions = position_matrix[atom_ids, :]
    centroid = get_centroid(position_matrix, atom_ids)
    return np.linalg.svd(atomic_positions - centroid)[-1][2, :]
