import pytest
from pytest_lazy_fixtures import lf

from .metal_topologies import *

# All fixtures must be visible for lazy_fixture() call.
from .three_plus_four import *
from .three_plus_three import *
from .two_plus_five import *
from .two_plus_four import *
from .two_plus_three import *


@pytest.fixture(
    params=(
        lf("cage_six_plus_eight"),
        lf("cage_four_plus_four"),
        lf("cage_one_plus_one"),
        lf("cage_two_plus_two"),
        lf("cage_twelve_plus_thirty"),
        lf("cage_eight_plus_sixteen"),
        lf("cage_five_plus_ten"),
        lf("cage_four_plus_eight"),
        lf("cage_six_plus_twelve"),
        lf("cage_ten_plus_twenty"),
        lf("cage_three_plus_six"),
        lf("cage_two_plus_four"),
        lf("cage_eight_plus_twelve"),
        lf("cage_four_plus_six"),
        lf("cage_four_plus_six_2"),
        lf("cage_six_plus_nine"),
        lf("cage_twenty_plus_thirty"),
        lf("cage_two_plus_three"),
        lf("metal_cage_m2l4_lantern"),
        lf("metal_cage_m3l3_triangle"),
        lf("metal_cage_m3l6"),
        lf("metal_cage_m4l4_square"),
        lf("metal_cage_m4l4_tetrahedron"),
        lf("metal_cage_m4l6_tetrahedron"),
        lf("metal_cage_m4l6_tetrahedron_spacer"),
        lf("metal_cage_m4l8"),
        lf("metal_cage_m4l8x"),
        lf("metal_cage_m4l82"),
        lf("metal_cage_m9l18"),
        lf("metal_cage_m6l2l3_prism"),
        lf("metal_cage_m6l12_cube"),
        lf("metal_cage_m8l6_cube"),
        lf("metal_cage_m12l24"),
        lf("metal_cage_m24l48"),
    ),
)
def cage(request):
    return request.param
