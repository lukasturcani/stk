import pytest
from pytest_lazy_fixtures import lf

from .hexagonal import *

# Fixtures need to be visible for lazy_fixture() calls.
from .honeycomb import *
from .kagome import *
from .linkerless_honeycomb import *
from .periodic_hexagonal import *
from .periodic_honeycomb import *
from .periodic_kagome import *
from .periodic_linkerless_honeycomb import *
from .periodic_square import *
from .square import *


@pytest.fixture(
    params=(
        lf("cof_honeycomb"),
        lf("cof_kagome"),
        lf("cof_square"),
        lf("cof_hexagonal"),
        lf("cof_linkerless_honeycomb"),
        lf("cof_periodic_honeycomb"),
        lf("cof_periodic_kagome"),
        lf("cof_periodic_square"),
        lf("cof_periodic_hexagonal"),
        lf("cof_periodic_linkerless_honeycomb"),
    ),
)
def cof(request):
    return request.param
