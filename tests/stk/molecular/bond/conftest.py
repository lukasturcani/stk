import pytest
import stk


@pytest.fixture(
    params=[
        stk.Bond(stk.H(0), stk.H(1), 1),
        stk.Bond(stk.He(2), stk.He(3), 2, (1, 0, -1)),
    ],
)
def bond(request):
    return request.param.clone()


@pytest.fixture(
    params=[
        lambda bond: None,
        lambda bond: {},
        lambda bond: {bond.atom1.id: stk.Be(33)},
        lambda bond: {bond.atom2.id: stk.Bi(122)},
        lambda bond: {
            bond.atom1.id: stk.K(4),
            bond.atom2.id: stk.S(7),
        },
    ],
)
def get_atom_map(request):
    return request.param
