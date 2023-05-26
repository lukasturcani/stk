import pytest
import stk

from ..case_data import CaseData

bb1 = stk.BuildingBlock("BrCCBr", [stk.BromoFactory()])
bb2 = stk.BuildingBlock("BrCC(CBr)CBr", [stk.BromoFactory()])
graph1 = stk.cage.FourPlusSix((bb1, bb2))

bb3 = stk.BuildingBlock("BrCNCBr", [stk.BromoFactory()])
bb4 = stk.BuildingBlock("BrCC(CNCBr)CBr", [stk.BromoFactory()])
graph2 = stk.cage.EightPlusTwelve((bb3, bb4))


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            crosser=stk.GeneticRecombination(
                get_gene=stk.BuildingBlock.get_num_functional_groups,
            ),
            records=(
                stk.MoleculeRecord(graph1),
                stk.MoleculeRecord(graph2),
            ),
            crossover_records=(
                stk.CrossoverRecord(
                    molecule_record=stk.MoleculeRecord(
                        topology_graph=graph1,
                    ),
                    crosser_name="GeneticRecombination",
                ),
                stk.CrossoverRecord(
                    molecule_record=stk.MoleculeRecord(
                        topology_graph=graph1.with_building_blocks(
                            building_block_map={bb1: bb3},
                        ),
                    ),
                    crosser_name="GeneticRecombination",
                ),
                stk.CrossoverRecord(
                    molecule_record=stk.MoleculeRecord(
                        topology_graph=graph1.with_building_blocks(
                            building_block_map={bb2: bb4},
                        ),
                    ),
                    crosser_name="GeneticRecombination",
                ),
                stk.CrossoverRecord(
                    molecule_record=stk.MoleculeRecord(
                        topology_graph=graph1.with_building_blocks(
                            building_block_map={
                                bb1: bb3,
                                bb2: bb4,
                            },
                        )
                    ),
                    crosser_name="GeneticRecombination",
                ),
                stk.CrossoverRecord(
                    molecule_record=stk.MoleculeRecord(
                        topology_graph=graph2.with_building_blocks(
                            building_block_map={
                                bb3: bb1,
                                bb4: bb2,
                            },
                        )
                    ),
                    crosser_name="GeneticRecombination",
                ),
                stk.CrossoverRecord(
                    molecule_record=stk.MoleculeRecord(
                        topology_graph=graph2.with_building_blocks(
                            building_block_map={bb4: bb2},
                        )
                    ),
                    crosser_name="GeneticRecombination",
                ),
                stk.CrossoverRecord(
                    molecule_record=stk.MoleculeRecord(
                        topology_graph=graph2.with_building_blocks(
                            building_block_map={bb3: bb1},
                        )
                    ),
                    crosser_name="GeneticRecombination",
                ),
                stk.CrossoverRecord(
                    molecule_record=stk.MoleculeRecord(
                        topology_graph=graph2,
                    ),
                    crosser_name="GeneticRecombination",
                ),
            ),
        ),
    ),
)
def genetic_recombination(request) -> CaseData:
    return request.param()
