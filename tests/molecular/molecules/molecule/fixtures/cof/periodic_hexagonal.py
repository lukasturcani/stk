import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicHexagonal(
                    building_blocks={
                        stk.BuildingBlock(
                            smiles='BrC1=C(Br)[C+]=N1',
                            functional_groups=[stk.BromoFactory()],
                        ): (
                            4, 5, 6, 7, 8, 9, 20, 21, 23, 24, 30, 36,
                            38, 40, 41, 42, 43, 46, 47, 52, 53, 60, 61,
                        ),
                        stk.BuildingBlock(
                            smiles='BrN1N(Br)[C+]=N1',
                            functional_groups=[stk.BromoFactory()],
                        ): (
                            10, 11, 12, 13, 14, 15, 22, 25, 26, 27, 28,
                            29, 37, 39, 44, 45, 54, 55, 56, 57, 58, 59,
                            31, 62, 63,
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                'Br[C+]1[C+]2[N+][C+2]C2(Br)[C+](I)[C+'
                                '](I)[C+](Br)[C+]1Br'
                            ),
                            functional_groups=[
                                stk.BromoFactory(),
                                stk.IodoFactory(),
                                stk.FluoroFactory(),
                            ],
                        ): (0, 1, 18, 50, 51),
                        stk.BuildingBlock(
                            smiles=(
                                'Br[C+]1[C+]2[S][C+2]C2(Br)[C+](I)[C+]'
                                '(I)[C+](Br)[C+]1Br'
                            ),
                            functional_groups=[
                                stk.BromoFactory(),
                                stk.IodoFactory(),
                                stk.FluoroFactory(),
                            ],
                        ): (2, 16, 34, 49),
                        stk.BuildingBlock(
                            smiles=(
                                'Br[C+]1[C+]2[S][O]C2(Br)[C+](I)[C+](I'
                                ')[C+](Br)[C+]1Br'
                            ),
                            functional_groups=[
                                stk.BromoFactory(),
                                stk.IodoFactory(),
                                stk.FluoroFactory(),
                            ],
                        ): (3, 17, 19, 32, 33, 35, 48),
                    },
                    lattice_size=(2, 2, 1),
                    vertex_alignments={0: 5},
                ),
            ),
            smiles=(
                '[C+]1=NC2=C1[C+]1[C+]3C4=C(N=[C+]4)[C+]4[C+]5[C+]6[C+'
                ']7[C+]8C9=C(N=[C+]9)[C+]9[C+]%10C%11=C([C+]=N%11)[C+]'
                '%11[C+]%12[C+]%13[C+]%14C%15=C(N=[C+]%15)[C+]%15[C+]%'
                '16C%17=C([C+]=N%17)[C+]%17[C+]([C+]%18[NH2+][C+2]C%18'
                '([C+]%18[C+]%19C%20=C(N=[C+]%20)C%20%21OS[C+]%20[C+]%'
                '20C%22=C(N=[C+]%22)[C+]%22[C+]([C+]%23C%24=C([C+]=N%2'
                '4)[C+]%24[C+]%25[C+]%26C%27=C([C+]=N%27)C%27%28OS[C+]'
                '%27[C+]%27C%29=C(N=[C+]%29)[C+]([C+]%29[C+]%20N%20[C+'
                ']=NN%20[C+]%20[C+]%30C%31=C(N=[C+]%31)[C+]%31[C+]2[C+'
                ']2C%32=C(N=[C+]%32)[C+]%32[C+]%33[C+](C%34=C(N=[C+]%3'
                '4)C4%34[C+2][NH2+][C+]8%34)[C+]4[C+]8[C+]%34SOC%34%32'
                'C%32=C(N=[C+]%32)[C+]%30[C+]%30C%32=C(N=[C+]%32)C%32%'
                '34OS[C+]%32[C+]%32C%35=C(N=[C+]%35)[C+]([C+]%16N%16N='
                '[C+]N%16[C+]%30[C+]%16S[C+2]C%16%20N%16N=[C+]N%16[C+]'
                '%17%19)[C+]%16SOC%16%17[C+]%15C%15=C([C+]=N%15)[C+]%1'
                '5[C+]%16C%19=C([C+]=N%19)C%19%20OS[C+]%19[C+]%19C%30='
                'C(N=[C+]%30)[C+]%30[C+]3N3[C+]=NN3[C+]([C+]%23N3N=[C+'
                ']N3[C+]%19[C+]3[C+](C%19=C(N=[C+]%19)[C+]%32[C+]%19[C'
                '+]([C+]%34N%23N=[C+]N8%23)N8[C+]=NN8[C+]([C+]9N8N=[C+'
                ']N48)[C+](N4N=[C+]N4[C+]%26[C+](N4N=[C+]N%194)C4([C+2'
                '][NH2+][C+]%244)N4N=[C+]N34)C3(OS[C+]%103)N3[C+]=NN3['
                'C+]%28[C+]([C+]3[C+]%27N4[C+]=NN4[C+]([C+]%31C4=C([C+'
                ']=N4)[C+]([C+]%16C4=C(N=[C+]4)C4([C+2][NH2+][C+]14)[C'
                '+]%30N1N=[C+]N%331)[C+]1S[C+2]C1([C+]%15N1N=[C+]N%131'
                ')N1N=[C+]N31)C1([C+2][NH2+][C+]21)N1N=[C+]N%291)N1[C+'
                ']=NN%121)[C+]%20N1N=[C+]N1%17)C1([C+2]S[C+]%221)N1[C+'
                ']=NN%181)[C+]%21N1N=[C+]N%251)N1N=[C+]N51)N1[C+]=NN61'
                ')N1[C+]=NN1[C+]%14C1([C+2]S[C+]%111)N1N=[C+]N71'
            ),
        ),
    ),
)
def cof_periodic_hexagonal(request):
    return request.param
