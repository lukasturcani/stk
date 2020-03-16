import stk

from ..utilities import is_equivalent_generic_functional_group


def test_repr(generic_functional_group):
    other = eval(repr(generic_functional_group), dict(stk.__dict__))
    is_equivalent_generic_functional_group(
        functional_group1=generic_functional_group,
        functional_group2=other,
    )
