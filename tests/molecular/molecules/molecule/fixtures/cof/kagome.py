import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Kagome(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrC1=C(Br)[C+]=N1',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+]('
                                'Br)[C+2]1'
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    lattice_size=(2, 2, 1),
                ),
            ),
            smiles=(
                'BrC1=C([C+]2[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5='
                'C([C+]=N5)[C+]5[C+2][C+]6C7=C(N=[C+]7)C37[C+]=N[C+]7['
                'C+]2C2=C([C+]=N2)[C+]2[C+]3[C+2][C+](C7=C(N=[C+]7)[C+'
                ']7[C+2][C+](C8=C(Br)[C+]=N8)[C+](C8=C(Br)[C+]=N8)[C+]'
                '8N=[C+]C78C7=C([C+]=N7)[C+]7[C+2][C+]8C9=C(N=[C+]9)C9'
                '%10[C+]=N[C+]9[C+](C9=C(N=[C+]9)[C+]9[C+]([C+2][C+]%1'
                '1C%12=C([C+]=N%12)[C+]%12[C+2][C+](C%13=C([C+]=N%13)C'
                '%13%14[C+]=N[C+]%13[C+](Br)[C+](Br)[C+2][C+]%14C%13=C'
                '(N=[C+]%13)C%11%13[C+]=N[C+]9%13)[C+](Br)[C+]9N=[C+]C'
                '%129C9=C([C+]=N9)[C+]5[C+]5N=[C+]C65Br)C5=C([C+]=N5)['
                'C+]8[C+]5N=[C+]C75C5=C3N=[C+]5)[C+]3[C+2][C+]%10C5=C('
                'N=[C+]5)[C+]5[C+2][C+](C6=C(Br)[C+]=N6)[C+](C6=C(Br)['
                'C+]=N6)[C+]6N=[C+]C56C5=C([C+]=N5)[C+]5[C+2][C+](Br)['
                'C+](Br)[C+]6N=[C+]C56C5=C3N=[C+]5)C3(C5=C(Br)[C+]=N5)'
                '[C+]=N[C+]23)[C+](Br)[C+]2N=[C+]C42C2=C(Br)[C+]=N2)N='
                '[C+]1'
            ),
        ),
    ),
)
def cof_kagome(request):
    return request.param
