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


@pytest.fixture(
    params=(
        CaseData(
            mutator=stk.RandomBuildingBlock(
                building_blocks=(bb2, ),
                is_replaceable=has_bromo,
            ),
            record=stk.ConstructedMoleculeRecord(
                molecule=stk.ConstructedMolecule(graph1),
                topology_graph=graph1,
            ),
            mutation_record=stk.ConstructedMoleculeMutationRecord(
                molecule_record=stk.ConstructedMoleculeRecord(
                    molecule=graph2,
                    topology_graph=graph2,
                ),
                mutator_name='RandomBuildingBlock',
            ),
        ),
    ),
)
def random_building_block(request):
    return request.param
