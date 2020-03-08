import pytest
from pytest_lazyfixture import lazy_fixture
import stk
import numpy as np

from .case_data import CaseData
# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=[
        stk.BuildingBlock.init(
            atoms=(stk.C(0, 4), ),
            bonds=(),
            position_matrix=np.array([[0.0, 0.0, 0.0]]),
            functional_groups=(),
        ),
        stk.BuildingBlock.init(
            atoms=(stk.C(0, 3), stk.H(1)),
            bonds=(stk.Bond(stk.C(0, 3), stk.H(1), 1), ),
            position_matrix=np.array([
                [0.39382080513175644, 0.0, 0.0],
                [-0.39382080513175644, 0.0, 0.0],
            ]),
            functional_groups=(),
        ),
        stk.BuildingBlock.init(
            atoms=(stk.C(0, 2), stk.H(1), stk.H(2)),
            bonds=(
                stk.Bond(stk.C(0, 2), stk.H(1), 1),
                stk.Bond(stk.C(0, 2), stk.H(2), 1),
            ),
            position_matrix=np.array([
                [-0.002271396061231665, 0.034037398527897535, -0.0],
                [-1.0494595365731274, -0.017073891221884126, -0.0],
                [1.0517309326343591, -0.016963507306017023, 0.0],
            ]),
            functional_groups=(),
        ),
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('N[C+][C+2]N'),
        stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()]),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='BrCNCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=2,
            ),
        ),
    ],
    scope='function',
)
def molecule(request):
    """
    A :class:`.Molecule` instance.

    """

    return request.param.clone()


def get_random_position_matrix(molecule):
    generator = np.random.RandomState(4)
    return generator.normal(
        loc=23.3,
        scale=32.1,
        size=(molecule.get_num_atoms(), 3),
    )


@pytest.fixture(
    params=(
        lambda molecule: np.zeros((molecule.get_num_atoms(), 3)),
        lambda molecule: np.array([
            [i, -i, 10.12*i] for i in range(molecule.get_num_atoms())
        ]),
        lambda molecule: molecule.get_position_matrix(),
        get_random_position_matrix,
    ),
)
def get_position_matrix(request):
    """
    Return a valid position matrix for a molecule.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which a position matrix is returned.

    Returns
    -------
    :class:`numpy.ndarray`
        A position matrix for `molecule`.

    """

    return request.param


@pytest.fixture(
    params=[
        [0, 0, 0],
        [10, 20, 30],
        [-10, 20, -30],
        [0.5, 10, -0.921],
    ]
)
def origin(request):
    return np.array(request.param)


@pytest.fixture
def get_origin_0(origin):
    return lambda molecule: origin


@pytest.fixture(
    params=[
        lambda molecule: molecule.get_centroid(),
        lazy_fixture('get_origin_0'),
    ],
)
def get_origin(request):
    """
    Return an origin parameter for a molecule.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule which needs an origin parameter for some method.

    Returns
    -------
    :class:`numpy.ndarray`
        The origin parameter to use.

    """

    return request.param


@pytest.fixture(
    params=[
        lambda molecule: None,
        lambda molecule: 0,
        lambda molecule: range(molecule.get_num_atoms()),
        lambda molecule: range(0, molecule.get_num_atoms(), 2),
        lambda molecule: list(
            range(0, min(1, molecule.get_num_atoms()))
        ),
        lambda molecule: tuple(
            range(0, min(1, molecule.get_num_atoms()))
        ),
        lambda molecule: (
            i for i in range(0, min(1, molecule.get_num_atoms()))
        ),
        lambda molecule: (),
        lambda molecule: range(min(molecule.get_num_atoms(), 1)),
        lambda molecule: range(min(molecule.get_num_atoms(), 2)),
        lambda molecule: range(min(molecule.get_num_atoms(), 3)),
    ],
)
def get_atom_ids(request):
    """
    Return an atom_ids parameter for a :class:`.Molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which `atom_ids` are returned.

    Retruns
    -------
    :class:`iterable` of :class:`int`
        An `atom_ids` parameter.

    """

    return request.param


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            smiles='[H]N([H])C([H])([H])C([H])([H])N([H])[H]',
        ),
        CaseData(
            molecule=stk.BuildingBlock('[H]NCCN'),
            smiles='[H]N([H])C([H])([H])C([H])([H])N([H])[H]',
        ),
        CaseData(
            molecule=stk.BuildingBlock('C(N)CN'),
            smiles='[H]N([H])C([H])([H])C([H])([H])N([H])[H]',
        ),
    ),
)
def building_block(request):
    return request.param


@pytest.fixture(
    params=(
        lazy_fixture('building_block'),
        lazy_fixture('cage'),
        lazy_fixture('cof'),
        lazy_fixture('polymer'),
        lazy_fixture('host_guest'),
        lazy_fixture('macrocycle'),
    ),
)
def case_data(request):
    return request.param
