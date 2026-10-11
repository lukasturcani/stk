import pytest
from pytest_lazy_fixtures import lf

# All fixtures must be visible for lazy_fixture() call.
from .octahedral import *
from .paddlewheel import *
from .porphyrin import *
from .square_planar import *


@pytest.fixture(
    params=(
        lf("metal_complex_octahedral"),
        lf("metal_complex_octahedral_lambda"),
        lf("metal_complex_octahedral_delta"),
        lf("metal_complex_porphyrin"),
        lf("metal_complex_paddlewheel"),
        lf("metal_complex_square_planar"),
        lf("metal_complex_bidentate_square_planar"),
        lf("metal_complex_cis_protected_square_planar"),
    ),
)
def metal_complex(request):
    return request.param
