from scipy.spatial.distance import euclidean


def get_maximum_diameter(position_matrix, atom_ids):
    atomic_positions = position_matrix[atom_ids, :]
    return float(
        euclidean(
            atomic_positions.min(axis=0),
            atomic_positions.max(axis=0),
        )
    )
