import pytest
import stk


@pytest.fixture(
    params=[
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


def test_with_atoms(bond, get_atom_map):
    bond.attr = 1
    bond._attr = 2
    atom_map = get_atom_map(bond)
    clone = bond.with_atoms(atom_map)
    assert clone is not bond

    expected_atom1 = atom_map.get(
        bond.get_atom1().get_id(),
        bond.get_atom1(),
    )
    assert clone.get_atom1() is expected_atom1

    expected_atom2 = atom_map.get(
        bond.get_atom2().get_id(),
        bond.get_atom2(),
    )
    assert clone.get_atom2() is expected_atom2

    assert bond.get_periodicity() == clone.get_periodicity()
    assert bond.get_order() == clone.get_order()
    assert bond.attr == clone.attr
    assert not hasattr(clone, '_attr')
