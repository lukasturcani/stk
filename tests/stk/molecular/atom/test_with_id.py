from .utilities import is_equivalent_atom


def test_with_id(atom, id):
    # Test immutability.
    before = atom.clone()
    _test_with_id(atom, id)
    is_equivalent_atom(before, atom)


def _test_with_id(atom, id):
    other = atom.with_id(id)
    assert other is not atom
    assert other.__class__ is atom.__class__
    assert other.get_id() == id
    assert other.get_charge() == atom.get_charge()
    assert other.get_atomic_number() == atom.get_atomic_number()
    assert other.get_mass() == atom.get_mass()
