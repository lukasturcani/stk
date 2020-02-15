import stk

from ..utilities import is_equivalent_functional_group


def test_repr(functional_group):
    other = eval(repr(functional_group), dict(stk.__dict__))
    is_equivalent_functional_group(functional_group, other)
