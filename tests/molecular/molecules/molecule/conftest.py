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
                    stk.BuildingBlock('BrC#CBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.polymer.Linear('AB', 2),
            ),
            smiles='BrC#C[C+]=NC#CC#C[C+]=NC#CBr',
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Honeycomb((2, 2, 1)),
            ),
            smiles=(
                'F[C+]1[C+](Br)[C+2][C+](C#CN=[C+][C+]2[C+2][C+]3[C+]'
                '=NC#C[C+]4[C+2][C+](C#CN=[C+][C+]5[C+2][C+]([C+]=NC#C'
                'Br)[C+](F)[C+]([C+]=NC#C[C+]6[C+2][C+](C#CN=[C+][C+]7'
                '[C+2][C+]([C+]=NC#CBr)[C+](F)[C+]([C+]=NC#CBr)[C+2]7)'
                '[C+2][C+](C#CN=[C+][C+]7[C+2][C+]([C+]=NC#C[C+]8[C+2]'
                '[C+](Br)[C+](F)[C+](C#CN=[C+][C+]([C+2]2)[C+]3F)[C+2]'
                '8)[C+2][C+]([C+]=NC#CBr)[C+]7F)[C+]6F)[C+2]5)[C+2][C+'
                '](Br)[C+]4F)[C+2][C+]1Br'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Honeycomb(
                    lattice_size=(2, 2, 1),
                    vertex_alignments={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
                ),
            ),
            smiles=(
                'F[C+]1[C+](Br)[C+2][C+](C#CN=[C+][C+]2[C+2][C+]3[C+]'
                '=NC#C[C+]4[C+2][C+](C#CN=[C+][C+]5[C+2][C+]([C+]=NC#C'
                'Br)[C+](F)[C+]([C+]=NC#C[C+]6[C+2][C+](C#CN=[C+][C+]7'
                '[C+2][C+]([C+]=NC#C[C+]8[C+2][C+](Br)[C+](F)[C+](C#CN'
                '=[C+][C+]([C+2]2)[C+]3F)[C+2]8)[C+2][C+]([C+]=NC#CBr)'
                '[C+]7F)[C+](F)[C+]([C+]=NC#C[C+]2[C+2][C+](C#CN=[C+]'
                'Br)[C+2][C+](C#CN=[C+]Br)[C+]2F)[C+2]6)[C+2]5)[C+2][C'
                '+](Br)[C+]4F)[C+2][C+]1Br'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Honeycomb(
                    lattice_size=(2, 2, 1),
                    vertex_alignments={0: 2, 1: 2},
                ),
            ),
            smiles=(
                'F[C+]1[C+](Br)[C+2][C+](C#CN=[C+][C+]2[C+2][C+]3[C+]'
                '=NC#C[C+]4[C+2][C+](C#CN=[C+][C+]5[C+2][C+]([C+]=NC#'
                'CBr)[C+](F)[C+]([C+]=NC#C[C+]6[C+2][C+](C#CN=[C+][C+'
                ']7[C+2][C+]([C+]=NC#CBr)[C+2][C+]([C+]=NC#CBr)[C+]7F)'
                '[C+](F)[C+](C#CN=[C+][C+]7[C+2][C+]([C+]=NC#C[C+]8[C'
                '+2][C+](Br)[C+](F)[C+](C#CN=[C+][C+]([C+2]2)[C+]3F)['
                'C+2]8)[C+2][C+]([C+]=NC#CBr)[C+]7F)[C+2]6)[C+2]5)[C+'
                '2][C+](Br)[C+]4F)[C+2][C+]1Br'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.LinkerlessHoneycomb((2, 2, 1)),
            ),
            smiles=(
                'F[C+]1[C+](Br)[C+2][C+]([C+]2[C+2][C+]3[C+](F)[C+]([C'
                '+2]2)[C+]2[C+2][C+]([C+2][C+](Br)[C+]2F)[C+]2[C+2][C'
                '+](Br)[C+](F)[C+]([C+2]2)[C+]2[C+2][C+]([C+]4[C+2][C+'
                '](Br)[C+](F)[C+](Br)[C+2]4)[C+2][C+]([C+]2F)[C+]2[C+2'
                '][C+]([C+2][C+](Br)[C+]2F)[C+]2[C+2][C+](Br)[C+](F)['
                'C+]3[C+2]2)[C+2][C+]1Br'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+](F)[C+](I)[C+](I)[C+](Br)'
                            'C1Br'
                        ),
                        functional_groups=[
                            stk.BromoFactory(),
                            stk.IodoFactory(),
                            stk.FluoroFactory(),
                        ],
                    ),
                ),
                topology_graph=stk.cof.Hexagonal(
                    lattice_size=(2, 2, 1),
                    vertex_alignments={0: 5},
                ),
            ),
            smiles=(
                'C1(Br)[C+]([C+]=NC#CBr)[C+]2[C+]=NC#C[C+]3[C+](I)[C+]'
                '([C+]=NC#CBr)[C+]4C#CN=[C+][C+]5[C+]([C+]=NC#C[C+]6[C'
                '+]7C#CN=[C+]C4[C+]3[C+]=NC#C[C+]3[C+]4C#CN=[C+][C+]2['
                'C+]2[C+]=NC#C[C+]8[C+]9C#CN=[C+]C%10[C+](C#CN=[C+][C+'
                ']21)[C+]([C+]=NC#CBr)[C+](I)[C+](I)[C+]%10[C+]=NC#C[C'
                '+]1[C+](F)[C+](I)[C+]2C#CN=[C+][C+]%10[C+](I)[C+](I)['
                'C+]%11[C+]=NC#C[C+]%12[C+](F)[C+](I)[C+](I)[C+]([C+]='
                'NC#CBr)C%12C#CN=[C+][C+]%12[C+](I)[C+]([C+]=NC#CBr)C%'
                '13[C+]=NC#C[C+]%14[C+](I)[C+]([C+]=NC#CBr)C%15C#CN=[C'
                '+][C+]%16[C+](I)[C+]([C+]=NC#CBr)C([C+]=NC#CBr)[C+]%1'
                '7C#CN=[C+][C+]%18[C+]([C+]=NC#CBr)[C+]([C+]=NC#CBr)[C'
                '+]%19[C+]=NC#C[C+]([C+]7[C+]=NC#CC3[C+]3[C+]=NC#CC%19'
                '[C+]%18[C+]=NC#C[C+]7[C+]%18[C+]=NC#C[C+]3[C+]4C#CN=['
                'C+]C8[C+]3[C+]=NC#C[C+]%18[C+]4C#CN=[C+][C+]8[C+]([C+'
                ']=NC#C[C+]%14[C+]%15C#CN=[C+][C+]4C7[C+]=NC#C[C+]%16%'
                '17)[C+]4[C+]=NC#C[C+]%13[C+]%12C#CN=[C+]C%11[C+]%10C#'
                'CN=[C+][C+]4C4C#CN=[C+][C+]2C1C#CN=[C+][C+]9[C+]3C#CN'
                '=[C+][C+]84)[C+]([C+]=NC#CBr)C6[C+]=NC#CBr)[C+]([C+]='
                'NC#CBr)[C+]([C+]=NC#CBr)C([C+]=NC#CBr)[C+]5Br'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+](Br)[C+](F)[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Kagome((2, 2, 1)),
            ),
            smiles=(
                'F[C+]1[C+](Br)[C+](Br)[C+2][C+]2C#CN=[C+][C+]3[C+]4[C'
                '+]=NC#C[C+]5[C+2][C+]([C+]=NC#C[C+]21)[C+](Br)[C+](F)'
                '[C+]5[C+]=NC#C[C+]1[C+2][C+]2C#CN=[C+][C+]5[C+]([C+]='
                'NC#C[C+]6[C+2][C+]([C+]=NC#C[C+]2[C+](F)[C+]1Br)[C+]('
                'Br)[C+](F)[C+]6[C+]=NC#CBr)[C+2][C+]([C+]=NC#CBr)[C+]'
                '([C+]=NC#C[C+]1[C+]2[C+]=NC#C[C+]6[C+](C#CN=[C+][C+]7'
                '[C+]([C+]=NC#C[C+]([C+2]2)[C+]([C+]=NC#CBr)[C+]1F)[C+'
                '2][C+]([C+]=NC#CBr)[C+]([C+]=NC#CBr)[C+]7F)[C+2][C+]1'
                'C#CN=[C+][C+]2[C+]7C#CN=[C+][C+]8[C+2][C+]([C+]=NC#CB'
                'r)[C+]([C+]=NC#CBr)[C+](F)[C+]8[C+]=NC#C[C+]8[C+2][C+'
                '](Br)[C+](Br)[C+](F)[C+]8C#CN=[C+][C+]([C+2]7)[C+](C#'
                'CN=[C+][C+]([C+]([C+]=NC#C[C+]1[C+]6F)[C+2]4)[C+]3F)['
                'C+]2F)[C+]5F'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)C(F)(Br)[C+]1Br',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Square((2, 2, 1)),
            ),
            smiles=(
                'FC1(Br)C(Br)=C2[C+]=NC#CC3(F)C(Br)=C([C+]=NC#CBr)[C+]'
                '3[C+]=NC#CC3=C([C+]=NC#CBr)[C+]([C+]=NC#CBr)C3(F)C#CN'
                '=[C+]C3=C(C#CN=[C+][C+]21)C(F)(Br)[C+]3[C+]=NC#CBr'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)C([C+2]F)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.FourPlusSix(),
            ),
            smiles='',
        ),
    ),
)
def case_data(request):
    return request.param
