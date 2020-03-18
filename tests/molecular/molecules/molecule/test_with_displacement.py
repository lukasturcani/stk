import numpy as np
import pytest


@pytest.fixture
def displacement(origin):
    return origin


def test_with_displacement(molecule, displacement):
    """
    Test :meth:`.Molecule.with_displacement`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The moleclue to test.

    displacement : :class:`numpy.ndarray`
        The displacement the clone should have.

    Returns
    -------
    None : :class:`NoneType`

    """

    new = molecule.with_displacement(displacement)
    assert np.allclose(
        a=molecule.get_position_matrix()+displacement,
        b=new.get_position_matrix(),
        atol=1e-32,
    )
