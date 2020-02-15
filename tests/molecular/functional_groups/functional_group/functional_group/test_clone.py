from .utilities import is_clone_functional_group


def test_clone(functional_group):
    clone = functional_group.clone()
    is_clone_functional_group(functional_group, clone)
