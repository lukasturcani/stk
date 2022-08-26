import pytest
from pytest_lazyfixture import lazy_fixture

from .metal_topologies import *  # noqa

# All fixtures must be visible for lazy_fixture() call.
from .three_plus_four import *  # noqa
from .three_plus_three import *  # noqa
from .two_plus_five import *  # noqa
from .two_plus_four import *  # noqa
from .two_plus_three import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("cage_six_plus_eight"),
        lazy_fixture("cage_four_plus_four"),
        lazy_fixture("cage_one_plus_one"),
        lazy_fixture("cage_two_plus_two"),
        lazy_fixture("cage_twelve_plus_thirty"),
        lazy_fixture("cage_eight_plus_sixteen"),
        lazy_fixture("cage_five_plus_ten"),
        lazy_fixture("cage_four_plus_eight"),
        lazy_fixture("cage_six_plus_twelve"),
        lazy_fixture("cage_ten_plus_twenty"),
        lazy_fixture("cage_three_plus_six"),
        lazy_fixture("cage_two_plus_four"),
        lazy_fixture("cage_eight_plus_twelve"),
        lazy_fixture("cage_four_plus_six"),
        lazy_fixture("cage_four_plus_six_2"),
        lazy_fixture("cage_six_plus_nine"),
        lazy_fixture("cage_twenty_plus_thirty"),
        lazy_fixture("cage_two_plus_three"),
        lazy_fixture("metal_cage_m2l4_lantern"),
        lazy_fixture("metal_cage_m3l3_triangle"),
        lazy_fixture("metal_cage_m3l6"),
        lazy_fixture("metal_cage_m4l4_square"),
        lazy_fixture("metal_cage_m4l4_tetrahedron"),
        lazy_fixture("metal_cage_m4l6_tetrahedron"),
        lazy_fixture("metal_cage_m4l6_tetrahedron_spacer"),
        lazy_fixture("metal_cage_m4l8"),
        lazy_fixture("metal_cage_m6l2l3_prism"),
        lazy_fixture("metal_cage_m6l12_cube"),
        lazy_fixture("metal_cage_m8l6_cube"),
        lazy_fixture("metal_cage_m12l24"),
        lazy_fixture("metal_cage_m24l48"),
    ),
)
def cage(request):
    return request.param
