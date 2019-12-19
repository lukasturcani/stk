import pytest
import stk

atoms = (stk.H(0), stk.F(1), stk.P(2), stk.H(3), stk.C(4), stk.K(5))


@pytest.fixture(
    params=[
        stk.FunctionalGroup(
            atoms=atoms,
            bonders=atoms[:2],
            deleters=atoms[-2:],
            fg_type=stk.functional_groups.fg_types['bromine'],
        ),
    ],
)
def functional_group(request):
    return request.param.clone()
