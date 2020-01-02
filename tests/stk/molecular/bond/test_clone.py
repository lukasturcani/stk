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
        lambda bond: {bond.get_atom1().get_id(): stk.Be(33)},
        lambda bond: {bond.get_atom2().get_id(): stk.Bi(122)},
        lambda bond: {
            bond.get_atom1().get_id(): stk.K(4),
            bond.get_atom2().get_id(): stk.S(7),
        },
    ],
)
def get_atom_map(request):
    """
    Return a valid `atom_map` parameter for a bond.

    Parameters
    ----------
    bond : :class:`.Bond`
        The bond for which the `atom_map` parameter needs to be
        created.

    Returns
    -------
    :class:`dict`
        A valid `atom_map` parameter for :meth:`.Bond.clone`.

    """

    return request.param


def test_clone(bond, get_atom_map):
    bond.attr = 1
    bond._attr = 2
    atom_map = get_atom_map(bond)
    clone = bond.clone(atom_map)

    if atom_map is None:
        atom_map = {}

    is_atom_clone(bond.get_atom1(), clone.get_atom1(), atom_map)
    is_atom_clone(bond.get_atom2(), clone.get_atom2(), atom_map)
    assert bond.get_periodicity() == clone.get_periodicity()
    assert bond.get_order() == clone.get_order()
    assert bond.attr == clone.attr
    assert not hasattr(clone, '_attr')


def is_atom_clone(atom, clone, atom_map):
    if atom.get_id() in atom_map:
        return is_atom_clone(clone, atom_map[atom.get_id()], {})
    else:
        assert atom is not clone
        assert atom.get_id() == clone.get_id()
        assert atom.get_charge() == clone.get_charge()
        assert atom.__class__ is clone.__class__
