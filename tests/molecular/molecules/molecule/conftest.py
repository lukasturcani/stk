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
                ''
            ),
        ),
    ),
)
def case_data(request):
    return request.param
