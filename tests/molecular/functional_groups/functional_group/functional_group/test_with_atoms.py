import pytest
import itertools as it
import stk

from ..utilities import is_clone_functional_group


def get_atom_map_0(functional_group):
    """
    Get an atom_map with a single atom.

    """

    # Make sure new_id is always valid, by making 1 larger than the
    # biggest on in the functional group. This prevents two atoms in
    # the functional group from having the same id.
    new_id = max(functional_group.get_atom_ids()) + 1
    atoms = (stk.Li(new_id), )
    return dict(zip(functional_group.get_atom_ids(), atoms))


@pytest.fixture(
    params=[
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
        :meth:`~.FunctionalGroup.with_atoms` on `functional_group`.

    """

    return request.param


def test_with_atoms(functional_group, get_atom_map):
    before = functional_group.clone()
    _test_with_atoms(functional_group, get_atom_map)
    # Test immutability.
    is_clone_functional_group(before, functional_group)


def _test_with_atoms(functional_group, get_atom_map):
    atom_map = get_atom_map(functional_group)
    clone = functional_group.with_atoms(atom_map)

    is_modified_sequence(
        atoms1=functional_group.get_atoms(),
        atoms2=clone.get_atoms(),
        atom_map=atom_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_atom_ids(),
        ids2=clone.get_atom_ids(),
        atom_map=atom_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_placer_ids(),
        ids2=clone.get_placer_ids(),
        atom_map=atom_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_core_atom_ids(),
        ids2=clone.get_core_atom_ids(),
        atom_map=atom_map,
    )


def is_modified_id_sequence(ids1, ids2, atom_map):
    for id1, id2 in it.zip_longest(ids1, ids2):
        if id1 in atom_map:
            assert id2 == atom_map[id1].get_id()
        else:
            assert id2 == id1


def is_modified_sequence(atoms1, atoms2, atom_map):
    for atom1, atom2 in it.zip_longest(atoms1, atoms2):
        assert atom2 is atom_map.get(atom1.get_id(), atom1)
