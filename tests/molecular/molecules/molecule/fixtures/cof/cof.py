import pytest
from pytest_lazyfixture import lazy_fixture

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
        lazy_fixture("cof_honeycomb"),
        lazy_fixture("cof_kagome"),
        lazy_fixture("cof_square"),
        lazy_fixture("cof_hexagonal"),
        lazy_fixture("cof_linkerless_honeycomb"),
        lazy_fixture("cof_periodic_honeycomb"),
        lazy_fixture("cof_periodic_kagome"),
        lazy_fixture("cof_periodic_square"),
        lazy_fixture("cof_periodic_hexagonal"),
        lazy_fixture("cof_periodic_linkerless_honeycomb"),
    ),
)
def cof(request):
    return request.param
