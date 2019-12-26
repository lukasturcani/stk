import itertools as it
import stk
import pytest


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
def test_get_atoms(molecule, bonds):
    for bond1, bond2 in it.zip_longest(molecule.get_bonds(), bonds):
        print(bond1, bond2)
        is_equivalent_bond(bond1, bond2)


def is_equivalent_bond(bond1, bond2):
    assert bond1.order == bond2.order
    is_equivalent_atom(bond1.atom1, bond2.atom1)
    is_equivalent_atom(bond1.atom2, bond2.atom2)


def is_equivalent_atom(atom1, atom2):
    assert atom1 is not atom2
    assert atom1.id == atom2.id
    assert atom1.charge == atom2.charge
    assert atom1.__class__ is atom2.__class__
