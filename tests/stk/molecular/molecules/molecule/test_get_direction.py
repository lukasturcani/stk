import numpy as np


def test_get_direction(molecule, get_atom_ids, direction):
    position_matrix = get_position_matrix(
        molecule=molecule,
        atom_ids=get_atom_ids(molecule),
        direction=direction,
    )
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.allclose(
        a=direction,
        b=molecule.get_direction(get_atom_ids(molecule)),
        atol=1e-32,
    )


def get_position_matrix(molecule, atom_ids, direction):
    """
    Create a position matrix with a specific `direction`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which the position matrix is created.

    atom_ids : :class:`iterable` of :class:`int`
        The ids of atoms which are to have a specific direction.
        If ``None``, then all atoms are used.

    direction : :class:`numpy.ndarray`
        The desired direction.

    Returns
    -------
    :class:`numpy.ndarray`
        A position matrix for `molecule`.

    """

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    elif not isinstance(atom_ids, (tuple, list)):
        atom_ids = tuple(atom_ids)

    position_matrix = molecule.get_position_matrix()
    for atom_id in atom_ids:
        position_matrix[atom_id] = direction*atom_id
    return position_matrix
