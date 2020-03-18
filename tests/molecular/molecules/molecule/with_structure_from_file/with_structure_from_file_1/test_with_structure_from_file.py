import numpy as np


def test_with_structure_from_file(
    molecule,
    get_position_matrix,
    path,
):
    """
    Test :meth:`.Molecule.with_structure_from_file`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    get_position_matrix : :class:`callable
        Takes a single parameter, `molecule`, and returns a valid
        position matrix for this test.

    path : :class:`str`
        A path into which the test structure is written.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save a copy of the position matrix, to ensure the original is
    # not modified by the test.
    position_matrix = molecule.get_position_matrix()
    _test_with_structure_from_file(
        molecule=molecule,
        get_position_matrix=get_position_matrix,
        path=path,
    )
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))


def _test_with_structure_from_file(
    molecule,
    get_position_matrix,
    path,
):
    """
    Test :meth:`.Molecule.with_structure_from_file`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    get_position_matrix : :class:`callable
        Takes a single parameter, `molecule`, and returns a valid
        position matrix for this test.

    path : :class:`str`
        A path into which the test structure is written.

    Returns
    -------
    None : :class:`NoneType`

    """

    position_matrix = get_position_matrix(molecule)
    molecule.with_position_matrix(position_matrix).write(path)
    loaded = molecule.with_structure_from_file(path)
    assert np.allclose(
        a=position_matrix,
        b=loaded.get_position_matrix(),
        atol=1e-4,
    )
