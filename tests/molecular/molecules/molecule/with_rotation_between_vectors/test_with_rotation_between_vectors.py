import numpy as np
import stk

from ...utilities import (
    get_displacement_vector,
    has_same_structure,
    is_clone,
)


def test_with_rotation_between_vectors(molecule, target, get_origin):
    """
    Test :meth:`.Molecule.with_rotation_between_vectors`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    target : :class:`numpy.ndarray`
        The target vector which defines the applied rotation.

    get_origin : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `origin` parameter for this test. Using a :class:`callable`
        allows the origin for the test to be set to molecule-specific
        values, such as the centroid of the molecule being tested.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Use to check that immutability is not violated.
    clone = molecule.clone()
    _test_with_rotation_between_vectors(molecule, target, get_origin)
    has_same_structure(molecule, clone)


def _test_with_rotation_between_vectors(molecule, target, get_origin):
    """
    Test :meth:`.Molecule.with_rotation_between_vectors`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    target : :class:`numpy.ndarray`
        The target vector which defines the applied rotation.

    get_origin : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `origin` parameter for this test. Using a :class:`callable`
        allows the origin for the test to be set to molecule-specific
        values, such as the centroid of the molecule being tested.

    Returns
    -------
    None : :class:`NoneType`

    """

    start = get_displacement_vector(molecule, 0, 1)
    new = molecule.with_rotation_between_vectors(
        start=start,
        target=target,
        origin=get_origin(molecule),
    )
    is_clone(new, molecule)
    result = get_displacement_vector(new, 0, 1)
    assert np.allclose(
        a=stk.normalize_vector(result),
        b=stk.normalize_vector(target),
        atol=1e-12,
    )
    assert abs(np.linalg.norm(start) - np.linalg.norm(result)) < 1e-14
