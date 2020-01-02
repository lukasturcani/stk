from .utilities import is_equivalent_atom


def test_clone(atom):
    atom.attr = 1
    atom._attr = 2
    clone = atom.clone()
    is_equivalent_atom(atom, clone)
    assert not hasattr(clone, '_attr')
    assert atom.attr == clone.attr
