import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Hexagonal(
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
                'BrC1=C([C+]2[C+]3C4=C([C+]=N4)[C+]4[C+](Br)[C+]5[C+]6'
                'C7=C([C+]=N7)[C+]7[C+]8[C+]9[C+]%10C%11=C(N=[C+]%11)['
                'C+]%11[C+]%12C%13=C(N=[C+]%13)[C+]%13[C+]([C+](C%14=C'
                '(Br)[C+]=N%14)[C+](C%14=C(Br)[C+]=N%14)[C+]%14C%15=C('
                '[C+]=N%15)[C+]%15[C+](C%16=C([C+]=N%16)C%16(OS[C+]7%1'
                '6)[C+]%10N7N=[C+]N7C%147OS[C+]%137)[C+](C7=C(N=[C+]7)'
                'C7([C+2][NH2+][C+]47)[C+]6N4N=[C+]N4Br)[C+](C4=C(N=[C'
                '+]4)[C+]2[C+](N2N=[C+]N2Br)C2(N4N=[C+]N4Br)[C+2][NH2+'
                '][C+]2[C+]3Br)[C+]2S[C+2]C2(N2N=[C+]N2Br)[C+]%15N2N=['
                'C+]N2Br)N2N=[C+]N2[C+]2[C+]3C4=C(N=[C+]4)C4(OS[C+]%12'
                '4)[C+]4[C+]6[C+]%11N7[C+]=NN7[C+]7[C+]%10[C+]%11C%12='
                'C([C+]=N%12)C%12%13OS[C+]%12[C+](C%12=C(N=[C+]%12)[C+'
                ']%12[C+](Br)[C+](Br)[C+]%14C%15=C(N=[C+]%15)[C+]%15[C'
                '+](N%16N=[C+]N%16Br)[C+](C%16=C([C+]=N%16)[C+]([C+]%1'
                '1N%11[C+]=NN%11[C+]%12C%11(C%12=C(Br)[C+]=N%12)OS[C+]'
                '%14%11)[C+]%11[NH2+][C+2]C%117N7N=[C+]N97)[C+]([C+](N'
                '7N=[C+]N57)C5(Br)[C+2]S[C+]%155)N5N=[C+]N85)[C+](Br)['
                'C+](Br)[C+]5[C+]%13N7N=[C+]N7C78OS[C+]7[C+]7C9=C([C+]'
                '=N9)[C+]9[C+]([C+](Br)[C+](I)[C+]%11N%12N=[C+]N%12[C+'
                ']%12[C+](Br)[C+](Br)[C+](I)[C+](N%13N=[C+]N%13Br)C%13'
                '([C+2][NH2+][C+]%12%13)N%12[C+]=NN%12[C+]%12[C+](I)[C'
                '+](C%13=C(Br)[C+]=N%13)C%13%14[C+2][NH2+][C+]%13[C+]('
                'C%13=C(N=[C+]%13)[C+]7[C+]([C+]([C+]8N7N=[C+]N%107)N7'
                'N=[C+]N67)N6N=[C+]N6[C+]6[C+](C7=C%14N=[C+]7)[C+](I)['
                'C+](C7=C(Br)[C+]=N7)C7(OS[C+]7[C+]6N6[C+]=NN46)C4=C(N'
                '=[C+]4)[C+]3[C+](I)[C+](N3N=[C+]N3Br)C3(N4N=[C+]N4Br)'
                '[C+2]S[C+]23)[C+]%12N2[C+]=NN2C%112[C+2]S[C+]92)N2N=['
                'C+]N52)N=[C+]1'
            ),
        ),
        # Non-planar linear BB.
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Hexagonal(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='C1(C(C(C(C(C1Br)Br)Br)Br)Br)Br',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='C1=C(Br)C=C2C=CC(Br)=CC2=C1CCCC',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    lattice_size=(2, 2, 1),
                ),
            ),
            smiles=(
                '[H]C1=C([H])C2=C([H])C(C3([H])C([H])(Br)C([H])(Br)C(['
                'H])(Br)C4([H])C5=C([H])C6=C(C([H])([H])C([H])([H])C(['
                'H])([H])C([H])([H])[H])C([H])=C(C([H])=C6C([H])=C5[H]'
                ')C5([H])C([H])(Br)C([H])(Br)C6([H])C7=C([H])C(C([H])('
                '[H])C([H])([H])C([H])([H])C([H])([H])[H])=C8C([H])=C('
                'C([H])=C([H])C8=C7[H])C7([H])C([H])(Br)C([H])(Br)C8(['
                'H])C9=C([H])C%10=C(C([H])([H])C([H])([H])C([H])([H])C'
                '([H])([H])[H])C([H])=C(C([H])=C%10C([H])=C9[H])C9([H]'
                ')C([H])(Br)C([H])(Br)C%10([H])C%11=C([H])C%12=C(C([H]'
                ')([H])C([H])([H])C([H])([H])C([H])([H])[H])C([H])=C(C'
                '([H])=C%12C([H])=C%11[H])C%11([H])C([H])(Br)C%12([H])'
                'C%13=C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])('
                '[H])[H])=C%14C([H])=C(C([H])=C([H])C%14=C%13[H])C%13('
                '[H])C([H])(Br)C%14([H])C%15=C([H])C%16=C(C([H])([H])C'
                '([H])([H])C([H])([H])C([H])([H])[H])C([H])=C(C([H])=C'
                '%16C([H])=C%15[H])C%15([H])C([H])(Br)C([H])(C%16=C([H'
                '])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])='
                'C%17C([H])=C(Br)C([H])=C([H])C%17=C%16[H])C([H])(C%16'
                '=C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])'
                '[H])=C%17C([H])=C(Br)C([H])=C([H])C%17=C%16[H])C([H])'
                '(C%16=C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])'
                '([H])[H])=C%17C([H])=C(Br)C([H])=C([H])C%17=C%16[H])C'
                '%15([H])C%15=C([H])C(C([H])([H])C([H])([H])C([H])([H]'
                ')C([H])([H])[H])=C%16C([H])=C(C([H])=C([H])C%16=C%15['
                'H])C%15([H])C%16([H])C%17=C([H])C%18=C(C([H])([H])C(['
                'H])([H])C([H])([H])C([H])([H])[H])C([H])=C(C([H])=C%1'
                '8C([H])=C%17[H])C%14([H])C([H])(C%14=C([H])C(C([H])(['
                'H])C([H])([H])C([H])([H])C([H])([H])[H])=C%17C([H])=C'
                '(Br)C([H])=C([H])C%17=C%14[H])C%13([H])C%13=C([H])C(C'
                '([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])=C%14C'
                '([H])=C(C([H])=C([H])C%14=C%13[H])C%13([H])C%14([H])C'
                '%17=C([H])C%18=C(C([H])([H])C([H])([H])C([H])([H])C(['
                'H])([H])[H])C([H])=C(C([H])=C%18C([H])=C%17[H])C%12(['
                'H])C([H])(C%12=C([H])C(C([H])([H])C([H])([H])C([H])(['
                'H])C([H])([H])[H])=C%17C([H])=C(C([H])=C([H])C%17=C%1'
                '2[H])C%12([H])C([H])(C%17=C([H])C%18=C(C([H])([H])C(['
                'H])([H])C([H])([H])C([H])([H])[H])C([H])=C(C([H])=C%1'
                '8C([H])=C%17[H])C9([H])C%10([H])C9=C([H])C(C([H])([H]'
                ')C([H])([H])C([H])([H])C([H])([H])[H])=C%10C([H])=C(B'
                'r)C([H])=C([H])C%10=C9[H])C9([H])C%10=C([H])C(C([H])('
                '[H])C([H])([H])C([H])([H])C([H])([H])[H])=C%17C([H])='
                'C(C([H])=C([H])C%17=C%10[H])C8([H])C7([H])C7=C([H])C('
                'C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])=C8C('
                '[H])=C(C([H])=C([H])C8=C7[H])C7([H])C8([H])C%10=C([H]'
                ')C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])=C'
                '%17C([H])=C(C([H])=C([H])C%17=C%10[H])C6([H])C5([H])C'
                '5=C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H]'
                ')[H])=C6C([H])=C(C([H])=C([H])C6=C5[H])C5([H])C([H])('
                'C6=C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H'
                '])[H])=C%10C([H])=C(C([H])=C([H])C%10=C6[H])C43[H])C('
                '[H])(Br)C([H])(C3=C([H])C(C([H])([H])C([H])([H])C([H]'
                ')([H])C([H])([H])[H])=C4C([H])=C(Br)C([H])=C([H])C4=C'
                '3[H])C3([H])C4=C([H])C(C([H])([H])C([H])([H])C([H])(['
                'H])C([H])([H])[H])=C6C([H])=C(C([H])=C([H])C6=C4[H])C'
                '4([H])C([H])(Br)C([H])(C6=C([H])C(C([H])([H])C([H])(['
                'H])C([H])([H])C([H])([H])[H])=C%10C([H])=C(Br)C([H])='
                'C([H])C%10=C6[H])C6([H])C%10=C([H])C%17=C(C([H])([H])'
                'C([H])([H])C([H])([H])C([H])([H])[H])C([H])=C(C([H])='
                'C%17C([H])=C%10[H])C%10([H])C([H])(Br)C([H])(C%17=C(['
                'H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])'
                '=C%18C([H])=C(Br)C([H])=C([H])C%18=C%17[H])C([H])(C%1'
                '7=C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H]'
                ')[H])=C%18C([H])=C(Br)C([H])=C([H])C%18=C%17[H])C%17('
                '[H])C%18=C([H])C%19=C(C([H])([H])C([H])([H])C([H])([H'
                '])C([H])([H])[H])C([H])=C(C([H])=C%19C([H])=C%18[H])C'
                '%18([H])C%19([H])C%20=C([H])C(C([H])([H])C([H])([H])C'
                '([H])([H])C([H])([H])[H])=C%21C([H])=C(C([H])=C([H])C'
                '%21=C%20[H])C%20([H])C%21([H])C%22=C([H])C(C([H])([H]'
                ')C([H])([H])C([H])([H])C([H])([H])[H])=C%23C([H])=C(C'
                '([H])=C([H])C%23=C%22[H])C([H])(C%14([H])C%14=C([H])C'
                '%22=C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]'
                ')C([H])=C(C([H])=C%22C([H])=C%14[H])C%12([H])C%12([H]'
                ')C%14=C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])'
                '([H])[H])=C%22C([H])=C(C([H])=C([H])C%22=C%14[H])C%21'
                '([H])C%14([H])C%21=C([H])C%22=C(C([H])([H])C([H])([H]'
                ')C([H])([H])C([H])([H])[H])C([H])=C(C([H])=C%22C([H])'
                '=C%21[H])C([H])(C7([H])C7=C([H])C(C([H])([H])C([H])(['
                'H])C([H])([H])C([H])([H])[H])=C%21C([H])=C(C([H])=C(['
                'H])C%21=C7[H])C9%12[H])C([H])(C7=C([H])C(C([H])([H])C'
                '([H])([H])C([H])([H])C([H])([H])[H])=C9C([H])=C(C([H]'
                ')=C([H])C9=C7[H])C4([H])C6([H])C4=C([H])C6=C(C([H])(['
                'H])C([H])([H])C([H])([H])C([H])([H])[H])C([H])=C(C([H'
                '])=C6C([H])=C4[H])C%14([H])C%20([H])C4=C([H])C(C([H])'
                '([H])C([H])([H])C([H])([H])C([H])([H])[H])=C6C([H])=C'
                '(C([H])=C([H])C6=C4[H])C%10%17[H])C8([H])C4=C([H])C(C'
                '([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])=C6C(['
                'H])=C(C([H])=C([H])C6=C4[H])C53[H])C3([H])C4=C([H])C('
                'C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])=C5C('
                '[H])=C(C([H])=C([H])C5=C4[H])C%19([H])C([H])(C4=C([H]'
                ')C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])=C'
                '5C([H])=C(C([H])=C([H])C5=C4[H])C([H])(C%16([H])C4=C('
                '[H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]'
                ')=C5C([H])=C(C([H])=C([H])C5=C4[H])C%133[H])C([H])(C3'
                '=C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])'
                '[H])=C4C([H])=C(Br)C([H])=C([H])C4=C3[H])C%15([H])C3='
                'C([H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])['
                'H])=C4C([H])=C(Br)C([H])=C([H])C4=C3[H])C([H])(C3=C(['
                'H])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])'
                '=C4C([H])=C(Br)C([H])=C([H])C4=C3[H])C%18([H])C3=C([H'
                '])C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])='
                'C4C([H])=C(Br)C([H])=C([H])C4=C3[H])C%11([H])C3=C([H]'
                ')C(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])=C'
                '4C([H])=C(Br)C([H])=C([H])C4=C3[H])=C([H])C(C([H])([H'
                '])C([H])([H])C([H])([H])C([H])([H])[H])=C2C([H])=C1Br'
            ),
        ),
        # One placer atom linear BB.
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Hexagonal(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='C1(C(C(C(C(C1Br)Br)Br)Br)Br)Br',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='C(Br)Br',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    lattice_size=(2, 2, 1),
                ),
            ),
            smiles=(
                '[H]C([H])(Br)C1([H])C([H])(Br)C2([H])C([H])([H])C3([H'
                '])C([H])(Br)C4([H])C([H])([H])C5([H])C([H])(Br)C6([H]'
                ')C([H])([H])C7([H])C([H])(Br)C([H])(Br)C8([H])C([H])('
                '[H])C9([H])C([H])(Br)C([H])(Br)C%10([H])C([H])([H])C%'
                '11([H])C([H])(Br)C([H])(Br)C%12([H])C([H])([H])C%13(['
                'H])C([H])(Br)C([H])(Br)C([H])(Br)C([H])(C([H])([H])Br'
                ')C%13([H])C([H])([H])C%13([H])C([H])(Br)C([H])(C([H])'
                '([H])Br)C%14([H])C([H])([H])C%15([H])C([H])(Br)C([H])'
                '(C([H])([H])Br)C%16([H])C([H])([H])C%17([H])C([H])(Br'
                ')C([H])(C([H])([H])Br)C([H])(C([H])([H])Br)C%18([H])C'
                '([H])([H])C%19([H])C([H])(C([H])([H])Br)C([H])(C([H])'
                '([H])Br)C%20([H])C([H])([H])C%21([H])C([H])(C([H])([H'
                '])Br)C([H])(C([H])([H])Br)C([H])(C([H])([H])C2([H])C('
                '[H])(C([H])([H])Br)C1([H])C([H])([H])Br)C1([H])C([H])'
                '([H])C3([H])C([H])(C([H])([H])Br)C4([H])C([H])([H])C2'
                '([H])C3([H])C([H])([H])C5([H])C([H])(C([H])([H])Br)C6'
                '([H])C([H])([H])C4([H])C([H])(C([H])([H])C8([H])C7([H'
                '])C([H])([H])Br)C5([H])C([H])([H])C9([H])C%10([H])C(['
                'H])([H])C6([H])C7([H])C([H])([H])C%11([H])C%12([H])C('
                '[H])([H])C%13([H])C%14([H])C([H])([H])C7([H])C7([H])C'
                '([H])([H])C%15([H])C%16([H])C([H])([H])C8([H])C([H])('
                'C([H])([H])C%17%18[H])C9([H])C([H])([H])C%19([H])C%20'
                '([H])C([H])([H])C([H])(C2([H])C([H])([H])C%211[H])C1('
                '[H])C([H])([H])C9([H])C2([H])C([H])([H])C([H])(C4([H]'
                ')C([H])([H])C31[H])C5([H])C([H])([H])C6([H])C7([H])C('
                '[H])([H])C82[H]'
            ),
        ),
        # Two placer atom linear BB.
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Hexagonal(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='C1(C(C(C(C(C1Br)Br)Br)Br)Br)Br',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    lattice_size=(2, 2, 1),
                ),
            ),
            smiles=(
                '[H]C([H])(Br)C([H])([H])C1([H])C([H])(Br)C([H])(Br)C('
                '[H])(Br)C2([H])C([H])([H])C([H])([H])C3([H])C([H])(Br'
                ')C([H])(Br)C4([H])C([H])([H])C([H])([H])C5([H])C([H])'
                '(Br)C([H])(Br)C6([H])C([H])([H])C([H])([H])C7([H])C(['
                'H])(Br)C([H])(Br)C8([H])C([H])([H])C([H])([H])C9([H])'
                'C([H])(Br)C%10([H])C([H])([H])C([H])([H])C%11([H])C(['
                'H])(Br)C%12([H])C([H])([H])C([H])([H])C%13([H])C([H])'
                '(Br)C([H])(C([H])([H])C([H])([H])Br)C([H])(C([H])([H]'
                ')C([H])([H])Br)C([H])(C([H])([H])C([H])([H])Br)C%13(['
                'H])C([H])([H])C([H])([H])C%13([H])C([H])(C([H])([H])C'
                '([H])([H])Br)C([H])(C([H])([H])C([H])([H])Br)C%14([H]'
                ')C([H])([H])C([H])([H])C%15([H])C([H])(C([H])([H])C(['
                'H])([H])Br)C([H])(C([H])([H])C([H])([H])Br)C%16([H])C'
                '([H])([H])C([H])([H])C%17([H])C([H])(C([H])([H])C([H]'
                ')([H])Br)C([H])(C([H])([H])C([H])([H])Br)C([H])(Br)C%'
                '18([H])C([H])([H])C([H])([H])C%19([H])C([H])(C([H])(['
                'H])C([H])([H])Br)C([H])(Br)C%20([H])C([H])([H])C([H])'
                '([H])C%21([H])C([H])(C([H])([H])C([H])([H])Br)C([H])('
                'Br)C([H])(C([H])([H])C([H])([H])C12[H])C1([H])C([H])('
                '[H])C([H])([H])C3([H])C4([H])C([H])([H])C([H])([H])C2'
                '([H])C3([H])C([H])([H])C([H])([H])C5([H])C6([H])C([H]'
                ')([H])C([H])([H])C4([H])C([H])(C([H])([H])C([H])([H])'
                'C7([H])C8([H])C([H])([H])C([H])([H])Br)C5([H])C([H])('
                '[H])C([H])([H])C9([H])C([H])(C([H])([H])C([H])([H])Br'
                ')C%10([H])C([H])([H])C([H])([H])C6([H])C7([H])C([H])('
                '[H])C([H])([H])C%11([H])C([H])(C([H])([H])C([H])([H])'
                'Br)C%12([H])C([H])([H])C([H])([H])C%13([H])C%14([H])C'
                '([H])([H])C([H])([H])C7([H])C7([H])C([H])([H])C([H])('
                '[H])C%15([H])C%16([H])C([H])([H])C([H])([H])C8([H])C('
                '[H])(C([H])([H])C([H])([H])C%18%17[H])C9([H])C([H])(['
                'H])C([H])([H])C%19([H])C%20([H])C([H])([H])C([H])([H]'
                ')C([H])(C2([H])C([H])([H])C([H])([H])C%211[H])C1([H])'
                'C([H])([H])C([H])([H])C9([H])C2([H])C([H])([H])C([H])'
                '([H])C([H])(C4([H])C([H])([H])C([H])([H])C31[H])C5([H'
                '])C([H])([H])C([H])([H])C6([H])C7([H])C([H])([H])C([H'
                '])([H])C82[H]'
            ),
        ),
    ),
)
def cof_hexagonal(request):
    return request.param
