import numpy as np

from ...utilities import is_clone


def test_with_centroid(molecule, get_atom_ids, centroid):
    """
    Test :meth:`.Molecule.with_centroid`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.with_centroid`. This
        allows the testing of different values of this parameter.

    centroid : :class:`numpy.ndarray`
        The desired centroid of the clone molecule.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save a copy of the position matrix, to ensure the original
    # molecule is not modified by the test, because it is meant to be
    # immutable.
    position_matrix = molecule.get_position_matrix()
    _test_with_centroid(molecule, get_atom_ids, centroid)
    assert np.all(
        np.equal(
            position_matrix,
            molecule.get_position_matrix(),
        )
    )


def _test_with_centroid(molecule, get_atom_ids, centroid):
    """
    Test :meth:`.Molecule.with_centroid`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.with_centroid`. This
        allows the testing of different values of this parameter.

    centroid : :class:`numpy.ndarray`
        The desired centroid of the clone molecule.

    Returns
    -------
    None : :class:`NoneType`

    """

    new = molecule.with_centroid(
        position=centroid,
        atom_ids=get_atom_ids(molecule),
    )
    is_clone(new, molecule)
    assert np.allclose(
        a=centroid,
        b=new.get_centroid(atom_ids=get_atom_ids(molecule)),
        atol=1e-8,
    )
