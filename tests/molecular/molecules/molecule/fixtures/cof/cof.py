import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .honeycomb import *  # noqa
from .kagome import *  # noqa
from .square import *  # noqa
from .hexagonal import *  # noqa
from .linkerless_honeycomb import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('cof_honeycomb'),
        lazy_fixture('cof_kagome'),
        lazy_fixture('cof_square'),
        lazy_fixture('cof_hexagonal'),
        lazy_fixture('cof_linkerless_honeycomb'),
    ),
)
def cof(request):
    return request.param
