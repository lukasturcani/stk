import numpy as np


def test_get_center_of_mass(
    molecule,
    get_atom_ids,
    get_position_matrix,
    center_of_mass,
):
    position_matrix = get_position_matrix(
        molecule=molecule,
        atom_ids=get_atom_ids(molecule),
        center_of_mass=center_of_mass,
    )
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.allclose(
        a=center_of_mass,
        b=molecule.get_center_of_mass(get_atom_ids(molecule)),
        atol=1e-32,
    )


def get_position_matrix(molecule, atom_ids, center_of_mass):
    """
    Create a position matrix with a specific `center_of_mass`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which the position matrix is created.

    atom_ids : :class:`iterable` of :class:`int`
        The ids of atoms which are to have a specific
        `center_of_mass`. If ``None``, then all atoms are used.

    center_of_mass : :class:`numpy.ndarray`
        The desired center of mass.

    Returns
    -------
    :class:`numpy.ndarray`
        A position matrix for `molecule`.

    """

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    elif not isinstance(atom_ids, (list, tuple)):
        atom_ids = tuple(atom_ids)

    position_matrix = molecule.get_position_matrix()
    # First, place all atoms on the center of mass.
    position_matrix[atom_ids, :] = center_of_mass
    # Displace pairs of atoms such that the center of mass
    # remains the same.
    masses = [a.mass for a in molecule.get_atoms()]
    generator = np.random.RandomState(4)
    for atom_id in range(0, molecule.get_num_atoms()-1, 2):
        displacement = generator.normal(scale=50, size=3)
        position_matrix[atom_id] += displacement / masses[atom_id]
        position_matrix[atom_id+1] -= displacement / masses[atom_id]
    return position_matrix
