import numpy as np


def test_write(molecule, get_position_matrix, path):
    """
    Test :meth:`.Molecule.write`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to write.

    get_position_matrix : :class:`callable
        Takes a single parameter, `molecule`, and returns a valid
        position matrix to use for this test.

    path : :class:`str`
        The path to which the test writes.

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
