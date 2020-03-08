import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData.init_constructed_molecule(
            topology_graph=stk.cof.Hexagonal(
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
                            'Br[C+]1[C+]2[S][C+2]C2(Br)[C+](I)[C+](I)'
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
                            'Br[C+]1[C+]2[S][O]C2(Br)[C+](I)[C+](I)[C'
                            '+](Br)[C+]1Br'
                        ),
                        functional_groups=[
                            stk.BromoFactory(),
                            stk.IodoFactory(),
                            stk.FluoroFactory(),
                        ],
                    ),
                ),
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
    ),
)
def cof_hexagonal(request):
    return request.param
