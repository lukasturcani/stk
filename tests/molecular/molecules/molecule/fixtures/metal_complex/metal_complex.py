import pytest
from pytest_lazyfixture import lazy_fixture

# All fixtures must be visible for lazy_fixture() call.
from .octahedral import *  # noqa
from .paddlewheel import *  # noqa
from .porphyrin import *  # noqa
from .square_planar import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("metal_complex_octahedral"),
        lazy_fixture("metal_complex_octahedral_lambda"),
        lazy_fixture("metal_complex_octahedral_delta"),
        lazy_fixture("metal_complex_porphyrin"),
        lazy_fixture("metal_complex_paddlewheel"),
        lazy_fixture("metal_complex_square_planar"),
        lazy_fixture("metal_complex_bidentate_square_planar"),
        lazy_fixture("metal_complex_cis_protected_square_planar"),
    ),
)
def metal_complex(request):
    return request.param
