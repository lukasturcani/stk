import pytest
from pytest_lazyfixture import lazy_fixture

# All fixtures must be visible for lazy_fixture() call.
from .octahedral import *  # noqa
from .square_planar import *  # noqa
from .porphyrin import *  # noqa
from .paddlewheel import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('metal_complex_octahedral'),
        lazy_fixture('metal_complex_octahedrallambda'),
        lazy_fixture('metal_complex_octahedraldelta'),
        lazy_fixture('metal_complex_porphyrin'),
        lazy_fixture('metal_complex_paddlewheel'),
        lazy_fixture('metal_complex_squareplanar'),
        lazy_fixture('metal_complex_bidentatesquareplanar'),
        lazy_fixture('metal_complex_cisprotectedsquareplanar'),
    ),
)
def metal_complex(request):
    return request.param
