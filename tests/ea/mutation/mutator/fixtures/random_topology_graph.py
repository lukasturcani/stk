import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            mutator=stk.RandomTopologyGraph(
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
