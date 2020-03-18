import numpy as np


def test_with_position_matrix(molecule, get_position_matrix):
    """
    Test :meth:`.Molecule.with_position_matrix`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    get_position_matrix : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        position matrix to use for this test.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save a copy of the position matrix, to ensure the original is
    # not modified by the test.
    position_matrix = molecule.get_position_matrix()
    _test_with_position_matrix(molecule, get_position_matrix)
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))


def _test_with_position_matrix(molecule, get_position_matrix):
    """
    Test :meth:`.Molecule.with_position_matrix`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    get_position_matrix : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        position matrix to use for this test.

    Returns
    -------
    None : :class:`NoneType`

    """

    position_matrix = get_position_matrix(molecule)
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))
