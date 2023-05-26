import pytest
import stk

from ..case_data import CaseData


def get_eight_plus_twelve(graph):
    return stk.cage.EightPlusTwelve(graph.get_building_blocks())


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            mutator=stk.RandomTopologyGraph(
                replacement_funcs=(get_eight_plus_twelve,),
            ),
            record=stk.MoleculeRecord(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="BrCC(CBr)CBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                ),
            ),
            mutation_record=stk.MutationRecord(
                molecule_record=stk.MoleculeRecord(
                    topology_graph=stk.cage.EightPlusTwelve(
                        building_blocks=(
                            stk.BuildingBlock(
                                smiles="BrCCBr",
                                functional_groups=[stk.BromoFactory()],
                            ),
                            stk.BuildingBlock(
                                smiles="BrCC(CBr)CBr",
                                functional_groups=[stk.BromoFactory()],
                            ),
                        ),
                    ),
                ),
                mutator_name="RandomTopologyGraph",
            ),
        ),
    ),
)
def random_topology_graph(request) -> CaseData:
    return request.param()
