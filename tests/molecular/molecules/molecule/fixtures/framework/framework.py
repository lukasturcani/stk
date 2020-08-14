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
        lazy_fixture('framework_honeycomb'),
        lazy_fixture('framework_kagome'),
        lazy_fixture('framework_square'),
        lazy_fixture('framework_hexagonal'),
        lazy_fixture('framework_linkerless_honeycomb'),
    ),
)
def framework(request):
    return request.param
