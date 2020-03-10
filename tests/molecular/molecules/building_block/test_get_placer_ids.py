import pytest
import stk
import itertools as it


class CaseData:
    def __init__(self, building_block, placer_ids):
        self.building_block = building_block
        self.placer_ids = placer_ids


@pytest.fixture(
    params=(
        CaseData(
            building_block=stk.BuildingBlock('Br[C+2][C+2]Br'),
            placer_ids=(0, 1, 2, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='BrCCCBr',
                functional_groups=[stk.BromoFactory()],
                placer_ids=(),
            ),
            placer_ids=(),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='BrCCCBr',
                functional_groups=[stk.BromoFactory()],
                placer_ids=(2, ),
            ),
            placer_ids=(2, ),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                functional_groups=[stk.IodoFactory()],
            ),
            placer_ids=(0, 1, 2, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                functional_groups=[stk.BromoFactory()],
            ),
            placer_ids=(1, 2),
        ),
    ),
)
def case_data(request):
    return request.param


def test_get_placer_ids(case_data):
    _test_get_placer_ids(
        building_block=case_data.building_block,
        placer_ids=case_data.placer_ids,
    )


def _test_get_placer_ids(building_block, placer_ids):
    for placer1, placer2 in it.zip_longest(
        building_block.get_placer_ids(),
        placer_ids,
    ):
        assert placer1 == placer2
