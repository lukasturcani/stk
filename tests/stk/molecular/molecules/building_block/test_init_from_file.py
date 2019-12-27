import pytest
import stk

from .utilities import is_equivalent_building_block


@pytest.fixture(
    params=[
        'building_block.mol',
        'building_block.pdb',
    ],
)
def filename(request):
    return request.param


def test_init_from_file(tmpdir, filename, building_block):
    path = str(tmpdir / filename)
    building_block.write(path)

    loaded = stk.BuildingBlock.init_from_file(
        path=path,
        functional_groups=building_block.get_functional_groups(),
    )
    is_equivalent_building_block(building_block, loaded)
