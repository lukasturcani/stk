import pytest
import stk
from .utilities import is_clone_functional_group


def get_atom_map_0(functional_group):
    """
    Get an atom_map with a single atom.

    """

    # Make sure new_id is always valid.
    new_id = max(functional_group.get_atom_ids()) + 1
    atoms = (stk.Li(new_id), )
    return dict(zip(functional_group.get_atom_ids(), atoms))


@pytest.fixture(
    params=[
        lambda functional_group: None,
        lambda functional_group: {},
        get_atom_map_0,
    ],
)
def get_atom_map(request):
    """
    Return a valid `atom_map` parameter for a functional group.

    Parameters
    ----------
    functional_group : :class:`.FunctionalGroup`
        A functional group being cloned, for which a valid
        `atom_map` parameter needs to created.

    Returns
    -------
    :class:`dict`
        A valid `atom_map` parameter for calling
        :meth:`~.FunctionalGroup.clone` on `functional_group`.

    """

    return request.param


def test_clone(functional_group, get_atom_map):
    functional_group.attr = 1
    functional_group._attr = 2
    atom_map = get_atom_map(functional_group)
    clone = functional_group.clone(atom_map)
    assert clone.attr == functional_group.attr
    assert not hasattr(clone, '_attr')
    is_clone_functional_group(functional_group, clone, atom_map)
