import pytest
from pytest_lazyfixture import lazy_fixture
import stk
import numpy as np

from .case_data import CaseData


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
            building_blocks=(
                stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                stk.BuildingBlock('BrCNCCBr', [stk.BromoFactory()]),
            ),
            topology_graph=stk.polymer.Linear('AB', 2),
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
            smiles='NCCN',
        ),
        CaseData(
            molecule=stk.BuildingBlock('[H]NCCN'),
            smiles='NCCN',
        ),
        CaseData(
            molecule=stk.BuildingBlock('C(N)CN'),
            smiles='NCCN',
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='BrCNCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.polymer.Linear('AB', 2),
            ),
            smiles='N(CCBr)CCCCCNCCCBr',
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCNCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='Brc1cc(Br)c(F)c(Br)c1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Honeycomb((2, 2, 1)),
            ),
            smiles=(
                'C1=C(Br)C(F)=C(Br)C=C1CCNCC1=CC2=C(F)C(=C1)CNCCC1=CC'
                '(=CC(Br)=C1F)CCNCC1=CC(=C(F)C(CNCCBr)=C1)CNCCC1=C(F)'
                'C(=CC(CCNCC3=CC(CNCCBr)=C(F)C(CNCCBr)=C3)=C1)CCNCC1=C'
                'C(=CC(CNCCBr)=C1F)CNCCC1=CC(=C(F)C(Br)=C1)CCNC2'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCNCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='Brc1cc(Br)c(F)c(Br)c1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Honeycomb(
                    lattice_size=(2, 2, 1),
                    vertex_alignments={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
                ),
            ),
            smiles=(
                'C1=C(Br)C(F)=C(Br)C=C1CCNCC1=CC2=C(F)C(=C1)CNCCC1=C'
                'C(=CC(Br)=C1F)CCNCC1=CC(=C(F)C(CNCCBr)=C1)CNCCC1=CC'
                '(=C(F)C(CNCCC3=CC(CCNCBr)=CC(CCNCBr)=C3F)=C1)CCNCC1='
                'CC(=CC(CNCCBr)=C1F)CNCCC1=CC(=C(F)C(Br)=C1)CCNC2'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCNCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='Brc1cc(Br)c(F)c(Br)c1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Honeycomb(
                    lattice_size=(2, 2, 1),
                    vertex_alignments={0: 2, 1: 2},
                ),
            ),
            smiles=(
                'C1=C(Br)C(F)=C(Br)C=C1CCNCC1=CC2=C(F)C(=C1)CNCCC1=C'
                'C(=CC(Br)=C1F)CCNCC1=CC(=C(F)C(CNCCBr)=C1)CNCCC1=C'
                'C(=CC(CCNCC3=C(F)C(CNCCBr)=CC(CNCCBr)=C3)=C1F)CCNC'
                'C1=CC(=CC(CNCCBr)=C1F)CNCCC1=CC(=C(F)C(Br)=C1)CCNC2'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Brc1cc(Br)c(F)c(Br)c1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.LinkerlessHoneycomb((2, 2, 1)),
            ),
            smiles=(
                'C1=C(Br)C(F)=C(Br)C=C1C1=CC2=C(F)C(=C1)C1=CC(=CC(Br)'
                '=C1F)C1=CC(=C(F)C(Br)=C1)C1=C(F)C(=CC(C3=CC(Br)=C(F)'
                'C(Br)=C3)=C1)C1=CC(=CC(Br)=C1F)C1=CC2=C(F)C(Br)=C1'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Brc1cc(Br)c(F)c(Br)c1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.LinkerlessHoneycomb(
                    lattice_size=(2, 2, 1),
                    vertex_alignments={0: 1, 1: 1},
                ),
            ),
            smiles=(
                'C1=C(Br)C=C(Br)C(F)=C1C1=C(F)C2=CC(=C1)C1=CC(=CC(Br)'
                '=C1F)C1=CC(=C(F)C(Br)=C1)C1=C(F)C(=CC(C3=CC(Br)=C(F)'
                'C(Br)=C3)=C1)C1=CC(=CC(Br)=C1F)C1=CC2=C(F)C(Br)=C1'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Brc1cc(Br)c(F)c(Br)c1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.LinkerlessHoneycomb(
                    lattice_size=(2, 2, 1),
                    vertex_alignments={0: 2, 1: 2},
                ),
            ),
            smiles=(
                'C1=C(Br)C=C(C2=CC3=CC(=C2F)C2=CC(=CC(Br)=C2F)C2=CC'
                '(=C(F)C(Br)=C2)C2=C(F)C(=CC(C4=CC(Br)=C(F)C(Br)=C4)='
                'C2)C2=CC(=CC(Br)=C2F)C2=CC3=C(F)C(Br)=C2)C(F)=C1Br'
            ),
        ),
    ),
)
def case_data(request):
    return request.param
