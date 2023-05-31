import numpy as np
import stk

from ...utilities import (
    get_displacement_vector,
    has_same_structure,
    is_clone,
)


def test_with_rotation_to_minimize_angle(molecule):
    """
    Test :meth:`.Molecule.with_rotation_to_minimize_angle`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Use to check that immutability is not violated.
    clone = molecule.clone()
    _test_with_rotation_to_minimize_angle(molecule)
    has_same_structure(molecule, clone)


def _test_with_rotation_to_minimize_angle(molecule):
    """
    Test :meth:`.Molecule.with_rotation_to_minimize_angle`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    start = get_displacement_vector(molecule, 0, 1)
    target = get_displacement_vector(molecule, 0, 2)
    new = molecule.with_rotation_to_minimize_angle(
        start=start,
        target=target,
        axis=stk.normalize_vector(np.cross(start, target)),
        origin=next(molecule.get_atomic_positions((0,))),
    )
    is_clone(new, molecule)

    result = get_displacement_vector(new, 0, 1)
    assert np.allclose(
        a=stk.normalize_vector(result),
        b=stk.normalize_vector(target),
        atol=1e-12,
    )
    assert abs(np.linalg.norm(start) - np.linalg.norm(result)) < 1e-15
