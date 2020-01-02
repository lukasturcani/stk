import pytest
import stk
import numpy as np


@pytest.fixture(
    params=[
        stk.BuildingBlock.init_from_dict({
            'class': 'BuildingBlock',
            'functional_groups': '[]',
            'position_matrix': [[0.0, 0.0, 0.0]],
            'atoms': '(C(0, charge=4),)',
            'bonds': [],
            'identity_key': '[C+4]',
        }),
        stk.BuildingBlock.init_from_dict({
            'class': 'BuildingBlock',
            'functional_groups': '[]',
            'position_matrix': [
                [0.39382080513175644, 0.0, 0.0],
                [-0.39382080513175644, 0.0, 0.0],
            ],
            'atoms': '(C(0, charge=3), H(1))',
            'bonds': [
                {
                    'atom1_id': 0,
                    'atom2_id': 1,
                    'order': 1,
                    'periodicity': [0, 0, 0],
                },
            ],
            'identity_key': '[H][C+3]'
        }),
        stk.BuildingBlock.init_from_dict({
            'class': 'BuildingBlock',
            'functional_groups': '[]',
            'position_matrix': [
                [-0.002271396061231665, 0.034037398527897535, -0.0],
                [-1.0494595365731274, -0.017073891221884126, -0.0],
                [1.0517309326343591, -0.016963507306017023, 0.0],
            ],
            'atoms': '(C(0, charge=2), H(1), H(2))',
            'bonds': [
                {
                    'atom1_id': 0,
                    'atom2_id': 1,
                    'order': 1,
                    'periodicity': [0, 0, 0],
                },
            ],
            'identity_key': '[H][C+2][H]'
        }),
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('N[C+][C+2]N'),
        stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()]),
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
        pytest.lazy_fixture('get_origin_0'),
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
