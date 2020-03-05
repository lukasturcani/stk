import stk
import pytest


from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='BrCNCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                topology_graph=stk.polymer.Linear('AB', 2),
            ),
            smiles='N(CCBr)CCCCCNCCCBr',
        ),
    ),
)
def case_data(request):
    return request.param
