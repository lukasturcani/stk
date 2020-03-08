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
                topology_graph=stk.cage.OnePlusOne(),
            ),
            smiles=(
                'F[C+]1[C+2]C23[C+2][C+]4[C+2][C+]([C+]12)C12[C+2][C+]'
                '(F)[C+]1[C+]3[C+2][C+]4[C+2]2'
            ),
        ),
    ),
)
def cage_one_plus_one(request):
    return request.param
