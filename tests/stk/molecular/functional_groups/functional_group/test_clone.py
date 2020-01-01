import pytest
import stk
import itertools as it
from .utilities import is_atom_clone


@pytest.fixture
def functional_group(test_case):
    """
    A :class:`.FunctionalGroup` instance.

    """

    return test_case.functional_group


def get_atom_map_0(functional_group):
    """
    Get an atom_map with a single atom.

    """

    atoms = (stk.Li(200), )
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
    is_clone_sequence(
        functional_group.get_atoms(),
        clone.get_atoms(),
        atom_map,
    )
    is_clone_sequence(
        functional_group.get_bonders(),
        clone.get_bonders(),
        atom_map,
    )
    is_clone_sequence(
        functional_group.get_deleters(),
        clone.get_deleters(),
        atom_map,
    )


def is_clone_sequence(atoms, clones, atom_map):
    """
    Test if atoms are correctly cloned.

    """

    if atom_map is None:
        atom_map = {}

    for atom, clone in it.zip_longest(atoms, clones):
        is_atom_clone(atom_map.get(atom.id, atom), clone)
