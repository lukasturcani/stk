import itertools as it
import stk
import pytest

from ..utilities import is_equivalent_bond


@pytest.mark.parametrize(
    argnames=('molecule', 'bonds'),
    argvalues=(
        (
            stk.BuildingBlock('NCCN'),
            (
                stk.Bond(stk.N(0), stk.C(1), 1),
                stk.Bond(stk.C(1), stk.C(2), 1),
                stk.Bond(stk.C(2), stk.N(3), 1),
                stk.Bond(stk.N(0), stk.H(4), 1),
                stk.Bond(stk.N(0), stk.H(5), 1),
                stk.Bond(stk.C(1), stk.H(6), 1),
                stk.Bond(stk.C(1), stk.H(7), 1),
                stk.Bond(stk.C(2), stk.H(8), 1),
                stk.Bond(stk.C(2), stk.H(9), 1),
                stk.Bond(stk.N(3), stk.H(10), 1),
                stk.Bond(stk.N(3), stk.H(11), 1),
            ),
        ),
    ),
)
def test_get_bonds(molecule, bonds):
    for bond1, bond2 in it.zip_longest(molecule.get_bonds(), bonds):
        is_equivalent_bond(bond1, bond2)
