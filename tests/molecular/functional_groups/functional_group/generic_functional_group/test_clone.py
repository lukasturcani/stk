from ..utilities import is_clone_generic_functional_group


def test_clone(generic_functional_group):
    clone = generic_functional_group.clone()
    is_clone_generic_functional_group(generic_functional_group, clone)
