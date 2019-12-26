import pytest
import numpy as np
import stk

from .utilities import has_same_structure, get_displacement_vector


@pytest.fixture(
    params=(
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('Brc1ccc(Br)cc1Br', [stk.BromoFactory()]),
    )
)
def molecule(request):
    """
    A :class:`.Molecule` instance which gets rotated.

    The molecule must have at least 2 atoms for the test to work.

    """

    return request.param.clone()


@pytest.fixture(
    params=[
        np.array([1., 0., 0.]),
        np.array([0., 1., 0.]),
        np.array([0., 0., 1.]),
        np.array([1., 1., 1.]),
    ],
)
def target(request):
    """
    The target vector onto which a molecule is rotated.

    """

    return request.param


def test_with_rotation_between_vectors(molecule, target, get_origin):
    # Use to check that immutability is not violated.
    clone = molecule.clone()
    _test_with_rotation_between_vectors(molecule, target, get_origin)
    has_same_structure(molecule, clone)


def _test_with_rotation_between_vectors(molecule, target, get_origin):
    start = get_displacement_vector(molecule, 0, 1)
    new = molecule.with_rotation_between_vectors(
        start=start,
        target=target,
        origin=get_origin(molecule),
    )
    result = get_displacement_vector(new, 0, 1)
    assert np.allclose(
        a=stk.normalize_vector(result),
        b=stk.normalize_vector(target),
        atol=1e-12,
    )
    assert abs(np.linalg.norm(start) - np.linalg.norm(result)) < 1e-32
