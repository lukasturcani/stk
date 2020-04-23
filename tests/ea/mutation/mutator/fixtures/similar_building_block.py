import pytest
import stk

from ..case_data import CaseData


def has_bromo(building_block):
    fg, = building_block.get_functional_groups(0)
    return fg.__class__ is stk.Bromo


bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
graph1 = stk.polymer.Linear((bb1, ), 'A', 2)

bb2 = stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()])
graph2 = stk.polymer.Linear((bb2, ), 'A', 2)

bb3 = stk.BuildingBlock('BrCNNCCNCBr', [stk.BromoFactory()])


@pytest.fixture(
    params=(
        CaseData(
            mutator=stk.SimilarBuildingBlock(
                building_blocks=(bb2, bb3),
                is_replaceable=has_bromo,
            ),
            record=stk.MoleculeRecord(graph1),
            mutation_record=stk.MutationRecord(
                molecule_record=stk.MoleculeRecord(graph2),
                mutator_name='SimilarBuildingBlock',
            ),
        ),
    ),
)
def similar_building_block(request):
    return request.param
