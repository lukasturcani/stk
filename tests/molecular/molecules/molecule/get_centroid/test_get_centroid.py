import numpy as np


def test_get_centroid(case_data, get_atom_ids):
    """
    Test :meth:`.Molecule.get_centroid`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the correct atom
        positions.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter. This allows the testing of different
        values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_centroid(
        molecule=case_data.molecule,
        centroid=get_centroid(
            position_matrix=case_data.position_matrix,
            atom_ids=get_atom_ids(case_data.molecule),
        ),
        get_atom_ids=get_atom_ids,
    )


def _test_get_centroid(molecule, centroid, get_atom_ids):
    """
    Test :meth:`.Molecule.get_centroid`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    centroid : :class:`numpy.ndarray`
        The correct centroid of `molecule`.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter. This allows the testing of different
        values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.allclose(
        a=centroid,
        b=molecule.get_centroid(get_atom_ids(molecule)),
        atol=1e-32,
    )


def get_centroid(position_matrix, atom_ids):
    if atom_ids is None:
        atom_ids = range(len(position_matrix))
    elif isinstance(atom_ids, int):
        atom_ids = (atom_ids, )
    else:
        atom_ids = len(atom_ids)

    return np.sum(position_matrix[atom_ids, :], axis=0) / len(atom_ids)
