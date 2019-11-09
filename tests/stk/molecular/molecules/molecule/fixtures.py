import pytest
import stk
import numpy as np


@pytest.fixture(params=[
    stk.BuildingBlock('NCCN'),
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


class GetAtomPositionsTestCase:
    """
    Represents a test case for :func:`test_get_atom_positions`.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule being tested.

    atom_ids : :class:`tuple` of :class:`int`
        The ids of atoms being tested.

    true_positions : :class:`np.ndarray`
        An ``[n, 3]`` matrix, where ``n`` is the number of atoms
        in :attr:`atom_ids`. For each atom in :attr:`atom_ids`
        this matrix holds the correct position.

    """

    __slots__ = ['molecule', 'atom_ids', 'true_positions']

    def __init__(self, molecule, atom_ids):
        self.molecule = molecule
        self.atom_ids = atom_ids
        self.true_positions = (
            molecule.get_position_matrix()[atom_ids, :]
        )


def get_atom_positions_test_case_1(molecule):
    return GetAtomPositionsTestCase(
        molecule=molecule,
        atom_ids=tuple(range(len(molecule.atoms))),
    )


def get_atom_positions_test_case_2(molecule):
    return GetAtomPositionsTestCase(
        molecule=molecule,
        atom_ids=tuple(range(0, len(molecule.atoms), 2)),
    )


def get_atom_positions_test_case_3(molecule):
    return GetAtomPositionsTestCase(
        molecule=molecule,
        atom_ids=(),
    )


@pytest.fixture(params=[
    get_atom_positions_test_case_1,
    get_atom_positions_test_case_2,
    get_atom_positions_test_case_3,
])
def get_atom_positions_test_case(request, molecule):
    return request.param(molecule)
