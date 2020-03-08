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
                topology_graph=stk.cage.TwoPlusTwo(),
            ),
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+]([C+]5[C+2][C+]6[C+]'
                '7[C+](F)[C+2]C7([C+2]5)[C+]5[C+2][C+]2[C+2]C2([C+2][C'
                '+](F)[C+]52)[C+]2[C+2][C+]4[C+2]C64[C+2][C+](F)[C+]24'
                ')[C+]13'
            ),
        ),
    ),
)
def cage_two_plus_two(request):
    return request.param
