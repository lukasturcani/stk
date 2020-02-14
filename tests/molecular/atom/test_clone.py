from .utilities import is_equivalent_atom


def test_clone(atom):
    clone = atom.clone()
    is_equivalent_atom(atom, clone)
