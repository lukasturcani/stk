import pytest
import stk
import numpy as np


@pytest.fixture(params=[
    stk.BuildingBlock('NCCN'),
    stk.BuildingBlock('NCCN').set_position_matrix(np.zeros((12, 3))),
])
def molecule(request):
    return request.param


@pytest.fixture(params=[
    [0, 0, 0],
    [10, 20, 30],
    [-10, 20, -30],
    [0.5, 10, -0.921],
])
def displacement(request):
    return request.param


@pytest.fixture(params=[
    -np.pi/2,
    np.pi/2,
])
def angle(request):
    return request.param


@pytest.fixture(params=[
    np.array([0, 1, 0]),
    np.array([1, 0, 0]),
    np.array([1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)])
])
def axis(request):
    return request.param


@pytest.fixture(params=[
    np.array([0, 0, 0]),
    np.array([-1.2, 10.12, 3]),
])
def origin(request):
    return request.param


# Values for get_atom_ids should never cause any dependent test
# function to fail.
never_fail_get_atom_ids = [
    lambda molecule: None,
    lambda molecule: range(len(molecule.atoms)),
    lambda molecule: range(0, len(molecule.atoms), 2),
    lambda molecule: range(0, min(1, len(molecule.atoms))),
    lambda molecule: list(range(0, min(1, len(molecule.atoms)))),
    lambda molecule: tuple(range(0, min(1, len(molecule.atoms)))),
    lambda molecule: (
        i for i in range(0, min(1, len(molecule.atoms)))
    ),
]


# Values for get_atom_ids which might cause some functions to fail.
# Functions which expect to fail depend on this fixture.
@pytest.fixture(params=[
    *never_fail_get_atom_ids,
    pytest.param(
        lambda molecule: (),
        marks=pytest.mark.xfail(strict=True, raises=ValueError),
    ),
])
def get_atom_ids(request):
    return request.param


# Values for get_atom_ids which might cause some functions to fail.
# Functions which do not expect to fail depend on this fixture.
@pytest.fixture(params=[
    *never_fail_get_atom_ids,
    lambda molecule: (),
])
def get_atom_ids_no_fail(request):
    return request.param
