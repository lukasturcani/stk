import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles=(
                            'Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]C2(B'
                            'r)[C+2]1'
                        ),
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.cage.FourPlusFour(),
            ),
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+]([C+]5[C+2][C+]6[C+]'
                '7[C+](F)[C+2]C7([C+2]5)[C+]5[C+2][C+]7[C+2]C8([C+2][C'
                '+](F)[C+]58)[C+]5[C+2][C+]8[C+2]C9([C+2][C+](F)[C+]59'
                ')[C+]5[C+2][C+]([C+2]C69[C+2][C+](F)[C+]59)C56[C+2][C'
                '+](F)[C+]5[C+]([C+2][C+]4[C+2]6)[C+]4[C+2][C+]([C+]5['
                'C+](F)[C+2]C85[C+2]4)C45[C+2][C+](F)[C+]4[C+]7[C+2][C'
                '+]2[C+2]5)[C+]13'
            ),
        ),
    ),
)
def cage_four_plus_four(request):
    return request.param
