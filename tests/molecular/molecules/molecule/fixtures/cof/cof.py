import pytest
from pytest_lazyfixture import lazy_fixture

from .hexagonal import *  # noqa

# Fixtures need to be visible for lazy_fixture() calls.
from .honeycomb import *  # noqa
from .kagome import *  # noqa
from .linkerless_honeycomb import *  # noqa
from .periodic_hexagonal import *  # noqa
from .periodic_honeycomb import *  # noqa
from .periodic_kagome import *  # noqa
from .periodic_linkerless_honeycomb import *  # noqa
from .periodic_square import *  # noqa
from .square import *  # noqa


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
