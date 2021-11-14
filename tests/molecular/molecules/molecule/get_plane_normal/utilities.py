import numpy as np

from ..utilities import get_centroid


def get_plane_normal(position_matrix, atom_ids):
    # In this case, Molecule.get_plane_normal should raise the error,
    # so this should should just fail silently, and the explicit
    # failure of Molecule.get_plane_normal will be caught.
    if len(atom_ids) == 0:
        return

    atomic_positions = position_matrix[atom_ids, :]
    centroid = get_centroid(position_matrix, atom_ids)
    return np.around(
        a=np.linalg.svd(atomic_positions - centroid)[-1][2, :],
        decimals=14,
    )
