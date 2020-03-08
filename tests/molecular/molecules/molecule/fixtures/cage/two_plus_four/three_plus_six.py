import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    params=(
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
    ),
)
def cage_three_plus_six(request):
    return request.param
