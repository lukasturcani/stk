import pytest
import stk


@pytest.fixture
def bond(atom1, atom2, order, periodicity):
    """
    A :class:`.Bond` instance.

    """

    return stk.Bond(atom1, atom2, order, periodicity)


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
    """
    A function which returns a valid `atom_map` parameter for a bond.

    """

    return request.param


def test_clone(bond, get_atom_map):
    bond.attr = 1
    bond._attr = 2
    atom_map = get_atom_map(bond)
    clone = bond.clone(atom_map)

    if atom_map is None:
        atom_map = {}

    is_atom_clone(bond.atom1, clone.atom1, atom_map)
    is_atom_clone(bond.atom2, clone.atom2, atom_map)
    assert bond.periodicity == clone.periodicity
    assert bond.order == clone.order
    assert bond.attr == clone.attr
    assert not hasattr(clone, '_attr')


def is_atom_clone(atom, clone, atom_map):
    if atom.id in atom_map:
        assert atom_map[atom.id] is clone
    else:
        assert atom is not clone
        assert atom.id == clone.id
        assert atom.charge == clone.charge
        assert atom.__class__ is clone.__class__
