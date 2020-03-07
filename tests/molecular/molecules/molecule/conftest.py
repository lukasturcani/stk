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
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='Br[C+]=NC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.polymer.Linear('AB', 2),
            ),
            smiles='BrC#CN=[C+]C1=C(C#CN=[C+]C2=C(Br)N=[C+]2)N=[C+]1',
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2[C+]=NC2(Br)[C+](F)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Honeycomb((2, 2, 1)),
            ),
            smiles=(
                '[H]C12[C+]=NC1(Br)[C+](F)[C+](Br)[C+2][C+]2C1=C([C+]2'
                '[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5=C(N=[C+]5)[C'
                '+]5[C+2][C+](C6=C(Br)[C+]=N6)[C+](F)C6(N=[C+]C56[H])C'
                '5=C([C+]=N5)C56N=[C+]C5([H])[C+](C5=C([C+]7[C+2][C+]('
                'C8=C(Br)[C+]=N8)[C+](F)C8(C9=C(Br)[C+]=N9)N=[C+]C78[H'
                '])N=[C+]5)[C+2][C+](C5=C(N=[C+]5)[C+]5[C+2][C+](C7=C('
                '[C+]=N7)[C+]7[C+2][C+](Br)[C+](F)C8(N=[C+]C78[H])C7=C'
                '(N=[C+]7)C7(N=[C+]C27[H])[C+]3F)C2([H])[C+]=NC2(C2=C('
                'Br)[C+]=N2)[C+]5F)[C+]6F)C2([H])[C+]=NC2(Br)[C+]4F)N='
                '[C+]1'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2[C+]=NC2(Br)[C+](F)[C+](Br)'
                            '[C+2]1'
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
                '[H]C12[C+]=NC1(Br)[C+](F)[C+](Br)[C+2][C+]2C1=C([C+]2'
                '[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5=C(N=[C+]5)[C'
                '+]5[C+2][C+](C6=C(Br)[C+]=N6)[C+](F)C6(N=[C+]C56[H])C'
                '5=C([C+]=N5)[C+]5[C+2][C+](C6=C([C+]7[C+2][C+](C8=C(B'
                'r)N=[C+]8)C8([H])[C+]=NC8(C8=C(Br)N=[C+]8)[C+]7F)[C+]'
                '=N6)[C+](F)C6(N=[C+]C56[H])C5=C(N=[C+]5)[C+]5[C+2][C+'
                '](C6=C([C+]=N6)[C+]6[C+2][C+](Br)[C+](F)C7(N=[C+]C67['
                'H])C6=C(N=[C+]6)C6(N=[C+]C26[H])[C+]3F)C2([H])[C+]=NC'
                '2(C2=C(Br)[C+]=N2)[C+]5F)C2([H])[C+]=NC2(Br)[C+]4F)N='
                '[C+]1'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2[C+]=NC2(Br)[C+](F)[C+](Br)'
                            '[C+2]1'
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
                '[H]C12[C+]=NC1(Br)[C+](F)[C+](Br)[C+2][C+]2C1=C([C+]2'
                '[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5=C(N=[C+]5)[C'
                '+]5[C+2][C+](C6=C(Br)[C+]=N6)[C+](F)C6(N=[C+]C56[H])C'
                '5=C([C+]=N5)[C+]5[C+2][C+](C6=C(N=[C+]6)[C+]6[C+2][C+'
                '](C7=C([C+]=N7)[C+]7[C+2][C+](Br)[C+](F)C8(N=[C+]C78['
                'H])C7=C(N=[C+]7)C7(N=[C+]C27[H])[C+]3F)C2([H])[C+]=NC'
                '2(C2=C(Br)[C+]=N2)[C+]6F)C2([H])[C+]=NC2(C2=C(C36N=[C'
                '+]C3([H])[C+](C3=C(Br)[C+]=N3)[C+2][C+](C3=C(Br)[C+]='
                'N3)[C+]6F)N=[C+]2)[C+]5F)C2([H])[C+]=NC2(Br)[C+]4F)N='
                '[C+]1'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2[C+]=NC2(Br)[C+](F)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.LinkerlessHoneycomb((2, 2, 1)),
            ),
            smiles=(
                '[H]C12[C+]=NC1(Br)[C+](F)[C+](Br)[C+2][C+]2[C+]1[C+2]'
                '[C+]2[C+]3[C+2][C+]([C+]4[C+2][C+](Br)[C+](F)C5(N=[C+'
                ']C45[H])C45N=[C+]C4([H])[C+]([C+]4[C+2][C+](Br)[C+](F'
                ')C6(Br)N=[C+]C46[H])[C+2][C+]([C+]4[C+2][C+]([C+]6[C+'
                '2][C+](Br)[C+](F)C7(N=[C+]C67[H])C6(N=[C+]C16[H])[C+]'
                '2F)C1([H])[C+]=NC1(Br)[C+]4F)[C+]5F)C1([H])[C+]=NC1(B'
                'r)[C+]3F'
            ),
        ),
        CaseData.init_constructed_molecule(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='BrC1=C(Br)[C+]=N1',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='BrN1N(Br)[C+]=N1',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+]2[N+][C+2]C2(Br)[C+](I)[C+](I)'
                        '[C+](Br)[C+]1Br'
                    ),
                    functional_groups=[
                        stk.BromoFactory(),
                        stk.IodoFactory(),
                        stk.FluoroFactory(),
                    ],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+]2[S][C+2]C2(Br)[C+](I)[C+](I)[C+]'
                        '(Br)[C+]1Br'
                    ),
                    functional_groups=[
                        stk.BromoFactory(),
                        stk.IodoFactory(),
                        stk.FluoroFactory(),
                    ],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+]2[S][O]C2(Br)[C+](I)[C+](I)[C+]'
                        '(Br)[C+]1Br'
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
            building_block_vertices={
                0: (
                    4, 5, 6, 7, 8, 9, 20, 21, 23, 24, 30, 36, 38,
                    40, 41, 42, 43, 46, 47, 52, 53, 60, 61
                ),
                1: (
                    10, 11, 12, 13, 14, 15, 22, 25, 26, 27, 28, 29,
                    37, 39, 44, 45, 54, 55, 56, 57, 58, 59, 31, 62, 63,
                ),
                2: (0, 1, 18, 50, 51),
                3: (2, 16, 34, 49),
                4: (3, 17, 19, 32, 33, 35, 48),
            },
            smiles=(
                'BrC1=C([C+]2[C+](Br)[C+](I)[C+]3C4=C([C+]=N4)[C+]4[C+'
                '](Br)[C+](Br)[C+]5[C+]6N7N=[C+]N7C78OS[C+]7[C+]7C9=C('
                '[C+]=N9)[C+]9[C+]([C+](Br)[C+](I)[C+]%10N%11N=[C+]N%1'
                '1[C+]%11[C+](Br)[C+](Br)[C+](I)[C+](N%12N=[C+]N%12Br)'
                'C%12([C+2][NH2+][C+]%11%12)N%11[C+]=NN%11[C+]%11[C+]('
                'I)[C+](C%12=C(Br)[C+]=N%12)C%12%13[C+2][NH2+][C+]%12['
                'C+](C%12=C(N=[C+]%12)[C+]7[C+]7[C+]%12[C+]8N8N=[C+]N8'
                '[C+]8[C+](C%14=C([C+]=N%14)C6%14OS[C+]4%14)[C+]4[C+]6'
                'C%14=C(N=[C+]%14)[C+]%14[C+](C%15=C([C+]=N%15)[C+]2[C'
                '+]2SOC23N2N=[C+]N42)[C+]2S[C+2]C2(Br)[C+](N2N=[C+]N2B'
                'r)[C+]2[C+]%14N3N=[C+]N3[C+]3[C+]4C%14=C(N=[C+]%14)[C'
                '+]%14[C+]([C+](Br)[C+](N%15N=[C+]N%15Br)[C+]%15C%16=C'
                '(N=[C+]%16)[C+]%16[C+](Br)[C+]%17[NH2+][C+2]C%17(C%17'
                '=C(Br)[C+]=N%17)[C+](N%17N=[C+]N%17Br)[C+](N%17N=[C+]'
                'N%17Br)[C+]%16C%16=C([C+]=N%16)[C+]%16[C+](C%17=C(N=['
                'C+]%17)C%14%17[C+2][NH2+][C+]%15%17)[C+]%14C%15=C([C+'
                ']=N%15)C%15(OS[C+]4%15)[C+]4[C+]%15C%17=C(N=[C+]%17)['
                'C+]%17[C+]%18C%19=C(N=[C+]%19)[C+]%19[C+]([C+](C%20=C'
                '(Br)[C+]=N%20)[C+](C%20=C(Br)[C+]=N%20)[C+](C%20=C([C'
                '+]=N%20)[C+]%14[C+](N%14N=[C+]N%14Br)C%14(N%20N=[C+]N'
                '%20Br)[C+2]S[C+]%16%14)C%14(OS[C+]%19%14)N%14[C+]=NN4'
                '%14)N4N=[C+]N4[C+]4[C+]%14C%16=C(N=[C+]%16)C%16(OS[C+'
                ']%18%16)[C+]([C+]([C+]%17N%16[C+]=NN%16[C+]8C8([C+2]['
                'NH2+][C+]68)N6N=[C+]N6[C+]%153)N3[C+]=NN%123)N3N=[C+]'
                'N3[C+]3[C+]([C+](C6=C%13N=[C+]6)[C+](I)[C+](C6=C(Br)['
                'C+]=N6)C6(OS[C+]36)C3=C(N=[C+]3)[C+]%14[C+](I)[C+](N3'
                'N=[C+]N3Br)C3(N6N=[C+]N6Br)[C+2]S[C+]43)N3[C+]=NN73)N'
                '3[C+]=NN23)[C+]%11N2[C+]=NN2C%102[C+2]S[C+]92)N2N=[C+'
                ']N52)N=[C+]1'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cof.Kagome((2, 2, 1)),
            ),
            smiles=(
                'BrC1=C([C+]2[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5='
                'C([C+]=N5)C56[C+]=N[C+]5[C+](Br)[C+](Br)[C+2][C+]6C5='
                'C(N=[C+]5)C35[C+]=N[C+]5[C+]2C2=C(Br)[C+]=N2)[C+]2C3='
                'C(N=[C+]3)[C+]3[C+]5[C+2][C+]6C7=C([C+]=N7)[C+]7[C+2]'
                '[C+](C8=C([C+]=N8)C89[C+]=N[C+]8[C+](Br)[C+](Br)[C+2]'
                '[C+]9C8=C(N=[C+]8)C68[C+]=N[C+]38)[C+](Br)[C+]3N=[C+]'
                'C73C3=C([C+]=N3)[C+]3[C+2][C+]6C7=C(N=[C+]7)C78[C+]=N'
                '[C+]7[C+](C7=C([C+]=N7)[C+]7[C+]9[C+2][C+](C%10=C(N=['
                'C+]%10)[C+]%10[C+2][C+](C%11=C(Br)[C+]=N%11)[C+](C%11'
                '=C(Br)[C+]=N%11)[C+]%11N=[C+]C%10%11C%10=C([C+]=N%10)'
                '[C+]%10[C+2][C+](C%11=C(N=[C+]%11)C4%11[C+]=N[C+]2%11'
                ')[C+](C2=C5N=[C+]2)[C+]2N=[C+]C%102C2=C9N=[C+]2)C2(C4'
                '=C(Br)[C+]=N4)[C+]=N[C+]72)[C+](C2=C(Br)[C+]=N2)[C+2]'
                '[C+]8C2=C([C+]=N2)[C+]2[C+2][C+](C4=C([C+]=N4)C64[C+]'
                '=N[C+]4[C+]3Br)[C+](Br)[C+]3N=[C+]C23C2=C(Br)[C+]=N2)'
                'N=[C+]1'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
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
                'FC1(Br)C2=C(C3=C([C+]=N3)C3(F)C(=C(C4=C(Br)[C+]=N4)[C'
                '+]3Br)C3=C([C+]=N3)[C+]3C(C4=C(Br)[C+]=N4)=C(C4=C(Br)'
                '[C+]=N4)C3(F)C3=C(N=[C+]3)C3=C(C4=C(Br)[C+]=N4)C(F)(B'
                'r)[C+]3C3=C2N=[C+]3)[C+]1Br'
            ),
        ),
        CaseData.init_constructed_molecule(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='BrC1=C(Br)[C+]=[C+]1',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='BrC1=C(Br)[C+]=N1',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+]2[C+][C+2]C2(Br)'
                        '[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]C2(Br)'
                        '[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            topology_graph=stk.cage.FourPlusSix(),
            building_block_vertices={
                0: (4, 5, 6, 7, 8),
                1: 9,
                2: (0, 1, 2),
                3: 3,
            },
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)[C+]5['
                'C+2][C+]6[C+2]C7([C+2][CH+][C+]57)C5=C([C+]=[C+]5)[C+'
                ']5[C+2][C+]([C+2]C7([C+2][CH+][C+]57)C5=C2[C+]=[C+]5)'
                'C2=C([C+]=[C+]2)[C+]2[C+2][C+](C5=C4[C+]=[C+]5)[C+]4['
                'CH+][C+2]C4([C+2]2)C2=C6[C+]=[C+]2)[C+]13'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]C2(B'
                            'r)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.TwoPlusThree(),
            ),
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)C56[C+'
                '2][C+]([C+2][C+](C7=C2[C+]=N7)[C+]5[C+](F)[C+2]6)C2=C'
                '4[C+]=N2)[C+]13'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]C2(B'
                            'r)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.FourPlusSix2(),
            ),
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)C56[C+'
                '2][C+]7[C+2][C+](C8=C(N=[C+]8)C89[C+2][C+]([C+2][C+]('
                'C%10=C([C+]=N%10)C%10%11[C+2][C+]([C+2][C+](C%12=C2[C'
                '+]=N%12)[C+]%10[C+](F)[C+2]%11)C2=C4[C+]=N2)[C+]8[C+]'
                '(F)[C+2]9)C2=C7[C+]=N2)[C+]5[C+](F)[C+2]6)[C+]13'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]C2(B'
                            'r)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.SixPlusNine(),
            ),
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)[C+]5['
                'C+2][C+]6[C+2]C7([C+2][C+](F)[C+]57)C5=C(N=[C+]5)[C+]'
                '5[C+2][C+]([C+2]C7([C+2][C+](F)[C+]57)C5=C4[C+]=N5)C4'
                '=C(N=[C+]4)[C+]4[C+2][C+]5C7=C([C+]=N7)[C+]7[C+2][C+]'
                '(C8=C2[C+]=N8)[C+]2[C+](F)[C+2]C2([C+2]7)C2=C(N=[C+]2'
                ')[C+]2[C+2][C+]([C+2]C7([C+2][C+](F)[C+]27)C2=C6[C+]='
                'N2)C2=C(N=[C+]2)C2([C+2]4)[C+2][C+](F)[C+]52)[C+]13'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]C2(B'
                            'r)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.EightPlusTwelve(),
            ),
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)C56[C+'
                '2][C+]7[C+2][C+](C8=C(N=[C+]8)C89[C+2][C+]%10[C+2][C+'
                '](C%11=C([C+]=N%11)C%11%12[C+2][C+]([C+2][C+](C%13=C2'
                '[C+]=N%13)[C+]%11[C+](F)[C+2]%12)C2=C(N=[C+]2)[C+]2[C'
                '+2][C+]%11C%12=C(N=[C+]%12)C%12%13[C+2][C+]([C+2][C+]'
                '(C%14=C([C+]=N%14)[C+]%14[C+2][C+](C%15=C7[C+]=N%15)['
                'C+]7[C+](F)[C+2]C7([C+2]%14)C7=C([C+]=N7)[C+]7[C+2][C'
                '+]([C+2]C%14([C+2][C+](F)[C+]7%14)C7=C4[C+]=N7)C4=C(N'
                '=[C+]4)C4([C+2]2)[C+2][C+](F)[C+]%114)[C+]%12[C+](F)['
                'C+2]%13)C2=C%10[C+]=N2)[C+]8[C+](F)[C+2]9)[C+]5[C+](F'
                ')[C+2]6)[C+]13'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]C2(B'
                            'r)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.TwentyPlusThirty(),
            ),
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)C56[C+'
                '2][C+]7[C+2][C+](C8=C([C+]=N8)C89[C+2][C+]%10[C+2][C+'
                '](C%11=C(N=[C+]%11)[C+]%11[C+2][C+]%12[C+2]C%13([C+2]'
                '[C+](F)[C+]%11%13)C%11=C([C+]=N%11)C%11%13[C+2][C+]%1'
                '4[C+2][C+](C%15=C(N=[C+]%15)C%15%16[C+2][C+]%17[C+2]['
                'C+](C%18=C([C+]=N%18)C%18%19[C+2][C+]([C+2][C+](C%20='
                'C(N=[C+]%20)[C+]%20[C+2][C+](C%21=C2[C+]=N%21)[C+]2[C'
                '+](F)[C+2]C2([C+2]%20)C2=C%12[C+]=N2)[C+]%18[C+](F)[C'
                '+2]%19)C2=C(N=[C+]2)C2%12[C+2][C+]%18[C+2][C+](C%19=C'
                '([C+]=N%19)[C+]%19[C+2][C+](C%20=C([C+]=N%20)[C+]%20['
                'C+2][C+]%21C%22=C(N=[C+]%22)[C+]%22[C+2][C+]%23[C+2]C'
                '%24([C+2][C+](F)[C+]%22%24)C%22=C([C+]=N%22)C%22%24[C'
                '+2][C+]([C+2][C+](C%25=C(N=[C+]%25)C%25%26[C+2][C+](['
                'C+2][C+](C%27=C([C+]=N%27)C%27%28[C+2][C+]%29[C+2][C+'
                '](C%30=C(N=[C+]%30)[C+]%30[C+2][C+]([C+2]C%31([C+2][C'
                '+](F)[C+]%30%31)C%30=C([C+]=N%30)C%30%31[C+2][C+]([C+'
                '2][C+](C%32=C(N=[C+]%32)C%32%33[C+2][C+]([C+2][C+](C%'
                '34=C([C+]=N%34)C%34([C+2]%20)[C+2][C+](F)[C+]%21%34)['
                'C+]%32[C+](F)[C+2]%33)C%20=C7N=[C+]%20)[C+]%30[C+](F)'
                '[C+2]%31)C7=C(N=[C+]7)C7%20[C+2][C+]([C+2][C+](C%21=C'
                '%10[C+]=N%21)[C+]7[C+](F)[C+2]%20)C7=C(N=[C+]7)C7%10['
                'C+2][C+]([C+2][C+](C%20=C%14[C+]=N%20)[C+]7[C+](F)[C+'
                '2]%10)C7=C%29[C+]=N7)C7=C%23N=[C+]7)[C+]%27[C+](F)[C+'
                '2]%28)[C+]%25[C+](F)[C+2]%26)C7=C%17N=[C+]7)[C+]%22[C'
                '+](F)[C+2]%24)C7=C%18N=[C+]7)[C+]7[C+](F)[C+2]C7([C+2'
                ']%19)C7=C4[C+]=N7)[C+]2[C+](F)[C+2]%12)[C+]%15[C+](F)'
                '[C+2]%16)[C+]%11[C+](F)[C+2]%13)[C+]8[C+](F)[C+2]9)[C'
                '+]5[C+](F)[C+2]6)[C+]13'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.TwoPlusFour(),
            ),
            smiles=(
                '[C+]1=NC2=C1C13[C+]=N[C+]1[C+]1C4=C(N=[C+]4)[C+]4[C+]'
                '2[C+2][C+]2C5=C([C+]=N5)[C+]3[C+2][C+]1C1=C(N=[C+]1)C'
                '21[C+]=N[C+]41'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.ThreePlusSix(),
            ),
            smiles=(
                '[C+]1=NC2=C1C13[C+]=N[C+]1[C+]1C4=C(N=[C+]4)[C+]4[C+]'
                '2[C+2][C+]2C5=C([C+]=N5)[C+]5[C+2][C+](C6=C(N=[C+]6)C'
                '26[C+]=N[C+]46)[C+]2C4=C([C+]=N4)[C+]1[C+2][C+]3C1=C('
                'N=[C+]1)C51[C+]=N[C+]21'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.FourPlusEight(),
            ),
            smiles=(
                '[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C(N=[C+]4)C45[C+]=N[C+]'
                '4[C+]4C6=C([C+]=N6)[C+]3[C+]3N=[C+]C13C1=C(N=[C+]1)[C'
                '+]1[C+]2[C+2][C+]2C3=C(N=[C+]3)[C+]3[C+2][C+](C6=C([C'
                '+]=N6)C26[C+]=N[C+]16)[C+]1C2=C([C+]=N2)[C+]4[C+2][C+'
                ']5C2=C(N=[C+]2)C32[C+]=N[C+]12'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.FivePlusTen(),
            ),
            smiles=(
                '[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C([C+]=N4)C45[C+]=N[C+]'
                '4[C+]4C6=C(N=[C+]6)[C+]3[C+]3N=[C+]C13C1=C(N=[C+]1)[C'
                '+]1[C+]2[C+2][C+]2C3=C(N=[C+]3)[C+]3[C+2][C+]6C7=C(N='
                '[C+]7)[C+]7[C+2][C+]8C9=C(N=[C+]9)[C+]5[C+2][C+]4C4=C'
                '([C+]=N4)C84[C+]=N[C+]4[C+]7C4=C([C+]=N4)C64[C+]=N[C+'
                ']4[C+]3C3=C([C+]=N3)C23[C+]=N[C+]13'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.SixPlusTwelve(),
            ),
            smiles=(
                '[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C(N=[C+]4)[C+]4[C+]5[C+'
                '2][C+]6C7=C(N=[C+]7)[C+]7[C+]8[C+2][C+]9C%10=C(N=[C+]'
                '%10)[C+]%10[C+2][C+]%11C%12=C(N=[C+]%12)[C+]%12[C+2]['
                'C+](C%13=C([C+]=N%13)C9%13[C+]=N[C+]7%13)[C+](C7=C([C'
                '+]=N7)C67[C+]=N[C+]47)[C+]4N=[C+]C%124C4=C([C+]=N4)C3'
                '4[C+]=N[C+]4[C+]1C1=C(N=[C+]1)[C+]%11[C+]1N=[C+]C%101'
                'C1=C(N=[C+]1)[C+]1[C+2][C+]2[C+](C2=C5[C+]=N2)[C+]2N='
                '[C+]C12C1=C8[C+]=N1'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.EightPlusSixteen(),
            ),
            smiles=(
                '[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C(N=[C+]4)[C+]4[C+2][C+'
                ']5C6=C([C+]=N6)[C+]6[C+2][C+]7C8=C([C+]=N8)C38[C+]=N['
                'C+]8[C+]1C1=C([C+]=N1)C13[C+]=N[C+]1[C+]1C8=C([C+]=N8'
                ')C89[C+]=N[C+]8[C+]8C%10=C([C+]=N%10)[C+]7[C+]7N=[C+]'
                'C67C6=C(N=[C+]6)C67[C+]=N[C+]6[C+]6C%10=C([C+]=N%10)['
                'C+]5[C+]5N=[C+]C45C4=C([C+]=N4)[C+]4[C+]2[C+2][C+]2C5'
                '=C([C+]=N5)[C+]3[C+2][C+]1C1=C(N=[C+]1)[C+]1[C+2][C+]'
                '(C3=C(N=[C+]3)C23[C+]=N[C+]43)[C+]2C3=C([C+]=N3)[C+]6'
                '[C+2][C+]7C3=C([C+]=N3)[C+]8[C+2][C+]9C3=C(N=[C+]3)C1'
                '3[C+]=N[C+]23'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)[C+]=N1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+](Br)'
                            '[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.TenPlusTwenty(),
            ),
            smiles=(
                '[C+]1=NC2=C1C13[C+]=N[C+]1[C+]1C4=C([C+]=N4)[C+]4[C+2'
                '][C+]5C6=C(N=[C+]6)C67[C+]=N[C+]6[C+]6C8=C(N=[C+]8)[C'
                '+]8[C+]9[C+2][C+]%10C%11=C(N=[C+]%11)[C+]%11[C+2][C+]'
                '%12C%13=C(N=[C+]%13)[C+]%13[C+2][C+](C%14=C([C+]=N%14'
                ')C%10%14[C+]=N[C+]8%14)[C+]8C%10=C([C+]=N%10)[C+]6[C+'
                '2][C+]7C6=C(N=[C+]6)[C+]6[C+]([C+2][C+]7C%10=C(N=[C+]'
                '%10)[C+]%12[C+]%10N=[C+]C%11%10C%10=C([C+]=N%10)[C+]%'
                '10[C+2][C+]%11C%12=C([C+]=N%12)[C+]%12[C+2][C+]%14C%1'
                '5=C(N=[C+]%15)C%15%16[C+]=N[C+]%15[C+](C%15=C([C+]=N%'
                '15)C4%15[C+]=N[C+]%15[C+]5C4=C(N=[C+]4)[C+]%14[C+]4N='
                '[C+]C%124C4=C9N=[C+]4)[C+]2[C+2][C+]%16C2=C([C+]=N2)C'
                '%112[C+]=N[C+]2[C+]%10C2=C([C+]=N2)[C+]3[C+2][C+]1C1='
                'C(N=[C+]1)C71[C+]=N[C+]61)C1=C(N=[C+]1)C%131[C+]=N[C+'
                ']81'
            ),
        ),
    ),
)
def case_data(request):
    return request.param
