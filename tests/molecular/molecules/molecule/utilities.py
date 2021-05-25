import numpy as np


def get_centroid(position_matrix, atom_ids):
    # In this case, Molecule.get_centroid should raise the error, so
    # this should should just fail silently, and the explicit failure
    # of Molecule.get_centroid will be caught.
    if len(atom_ids) == 0:
        return

    return np.sum(position_matrix[atom_ids, :], axis=0) / len(atom_ids)


def get_direction(position_matrix, atom_ids):
    # In this case, Molecule.get_direction should raise the error, so
    # this should should just fail silently, and the explicit failure
    # of Molecule.get_direction will be caught.
    if len(atom_ids) == 0:
        return

    atom_positions = position_matrix[atom_ids, :]
    centered_positions = atom_positions - atom_positions.mean(axis=0)
    return np.around(np.linalg.svd(centered_positions)[-1][0], 14)
