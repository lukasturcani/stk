import numpy as np
import stk


def test_get_maximum_diameter(
    molecule,
    get_atom_ids,
    maximum_diameter,
):
    position_matrix = get_position_matrix(
        molecule=molecule,
        atom_ids=get_atom_ids(molecule),
        maximum_diameter=maximum_diameter,
    )
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.allclose(
        a=maximum_diameter,
        b=molecule.get_maximum_diameter(get_atom_ids(molecule)),
        atol=1e-32,
    )


def get_position_matrix(molecule, atom_ids, maximum_diameter):
    """
    Create a position matrix with a specific `maximum_diameter`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which the position matrix is created.

    atom_ids : :class:`iterable` of :class:`int`
        The ids of atoms which are to have a specific maximum
        diameter. If ``None``, then all atoms are used.

    maximum_diameter : :class:`float`
        The desired maximum diameter.

    Returns
    -------
    :class:`numpy.ndarray`
        The position matrix for `molecule`.

    """

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    elif not isinstance(molecule, (list, tuple)):
        atom_ids = tuple(atom_ids)

    position_matrix = molecule.get_position_matrix()
    position_matrix[atom_ids, :] = 0
    direction = get_direction_vector()
    atom1_id, atom2_id = next(atom_ids), next(atom_ids)
    position_matrix[atom1_id] += 3*direction*maximum_diameter/5
    position_matrix[atom2_id] -= 2*direction*maximum_diameter/5


def get_direction_vector():
    generator = np.random.RandomState(4)
    return stk.normalize_vector(generator.rand(3))
