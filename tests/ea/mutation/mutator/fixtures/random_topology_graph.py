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
                topology_graph=stk.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCC(CBr)CBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                ),
            ),
            mutation_record=stk.ConstructedMoleculeMutationRecord(
                molecule_record=stk.ConstructedMoleculeRecord(
                    topology_graph=stk.EightPlusTwelve(
                        building_blocks=(
                            stk.BuildingBlock(
                                smiles='BrCCBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                            stk.BuildingBlock(
                                smiles='BrCC(CBr)CBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                        ),
                    ),
                ),
                mutator_name='RandomTopologyGraph',
            ),
        ),
    ),
)
def random_topology_graph(request):
    return request.param
