import pytest
import stk
from ..case_data import CaseData


def has_bromo(building_block):
    fg, = building_block.get_functional_groups(0)
    return type(fg) is stk.Bromo


def _get_case_data_1() -> CaseData:

    bb1 = stk.BuildingBlock('BrCC(Cl)CBr', [stk.BromoFactory()])
    graph1 = stk.polymer.Linear((bb1, ), 'A', 2)

    bb2 = stk.BuildingBlock('BrCC(F)CBr', [stk.BromoFactory()])
    graph2 = stk.polymer.Linear((bb2, ), 'A', 2)

    return CaseData(
        mutator=stk.SubstituteSubstructure(
            query_smarts='[Cl]',
            replacement_smiles='[F]',
            is_replaceable=has_bromo,
            replacement_functional_groups=[stk.BromoFactory()],
        ),
        record=stk.MoleculeRecord(graph1),
        mutation_record=stk.MutationRecord(
            molecule_record=stk.MoleculeRecord(graph2),
            mutator_name='SubstituteSubstructure',
        ),
    )


@pytest.fixture(
    scope='session',
    params=(
        _get_case_data_1,
    ),
)
def substitute_substructure(request) -> CaseData:
    return request.param()
