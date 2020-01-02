from .utilities import is_clone_functional_group


def test_clone(functional_group):
    functional_group.attr = 1
    functional_group._attr = 2
    clone = functional_group.clone()
    assert clone.attr == functional_group.attr
    assert not hasattr(clone, '_attr')
    is_clone_functional_group(functional_group, clone)
