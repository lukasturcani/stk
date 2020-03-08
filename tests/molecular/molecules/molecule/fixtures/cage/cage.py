import pytest
from pytest_lazyfixture import lazy_fixture

# All fixtures must be visible for lazy_fixture() call.
from .three_plus_four import *  # noqa
from .three_plus_three import *  # noqa
from .two_plus_five import *  # noqa
from .two_plus_four import *  # noqa
from .two_plus_three import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('cage_six_plus_eight'),

    ),
)
def cage(request):
    return request.param
