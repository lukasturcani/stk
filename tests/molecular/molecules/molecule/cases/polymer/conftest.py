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
            smiles=(
                '[H]N(C([H])([H])C([H])([H])Br)C([H])([H])C([H])([H])'
                'C([H])([H])C([H])([H])C([H])([H])N([H])C([H])([H])C'
                '([H])([H])C([H])([H])Br'
            ),
        ),
    ),
)
def case_data(request):
    return request.param
