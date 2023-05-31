import pytest
import stk

from ..case_data import CaseData


def has_bromo(building_block):
    (fg,) = building_block.get_functional_groups(0)
    return fg.__class__ is stk.Bromo


def _get_case_data_1() -> CaseData:
    bb1 = stk.BuildingBlock("BrCCBr", [stk.BromoFactory()])
    graph1 = stk.polymer.Linear((bb1,), "A", 2)

    bb2 = stk.BuildingBlock("BrCNCBr", [stk.BromoFactory()])
    graph2 = stk.polymer.Linear((bb2,), "A", 2)

    return CaseData(
        mutator=stk.RandomMutator(
            mutators=(
                stk.RandomBuildingBlock(
                    building_blocks=(bb2,),
                    is_replaceable=has_bromo,
                ),
            ),
        ),
        record=stk.MoleculeRecord(graph1),
        mutation_record=stk.MutationRecord(
            molecule_record=stk.MoleculeRecord(graph2),
            mutator_name="RandomBuildingBlock",
        ),
    )


@pytest.fixture(
    scope="session",
    params=(_get_case_data_1,),
)
def random_mutator(request) -> CaseData:
    return request.param()
