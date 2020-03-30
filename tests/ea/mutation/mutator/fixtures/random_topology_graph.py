import pytest
import stk

from ...case_data import CaseData


def get_eight_plus_twelve(graph):
    return stk.EightPlusTwelve(graph.get_building_blocks())


@pytest.fixture(
    params=(
        CaseData(
            mutator=stk.RandomTopologyGraph(
                replacement_funcs=(get_eight_plus_twelve, ),
            ),
            record=stk.ConstructedMoleculeRecord(
            ),
            mutation_record=stk.ConstructedMoleculeMutationRecord(
            ),
        ),
    ),
)
def random_topology_graph(request):
    return request.param
