import stk

from .utilities import is_equivalent_atom


def test_repr(atom):
    other = eval(repr(atom), dict(stk.__dict__))
    is_equivalent_atom(other, atom)
