import itertools as it
import stk

from typing import Callable

from ...utilities import is_clone_functional_group


def test_with_ids(
    functional_group: stk.FunctionalGroup,
    get_id_map: Callable[[stk.FunctionalGroup], dict[int, int]],
) -> None:
    """
    Test :meth:`.FunctionalGroup.with_ids`.

    Parameters:

        functional_group:
            The functional group to test.

        get_id_map:
            Takes a single parameter, `functional_group`, and returns a
            valid `id_map` parameter for its
            :meth:`.FunctionalGroup.with_atoms` method. This allows the
            testing of different values of this parameter.

    """

    # Save a clone to ensure that "functional_group" is not changed by
    # the test, because it should be immutable.
    before = functional_group.clone()
    _test_with_atoms(functional_group, get_id_map)
    is_clone_functional_group(before, functional_group)


def _test_with_atoms(
    functional_group: stk.FunctionalGroup,
    get_id_map: Callable[[stk.FunctionalGroup], dict[int, int]],
) -> None:
    """
    Test :meth:`.FunctionalGroup.with_atoms`.

    Parameters:

        functional_group:
            The functional group to test.

        get_id_map:
            Takes a single parameter, `functional_group`, and returns a
            valid `id_map` parameter for its
            :meth:`.FunctionalGroup.with_atoms` method. This allows the
            testing of different values of this parameter.

    """

    id_map = get_id_map(functional_group)
    clone = functional_group.with_ids(id_map)

    is_modified_id_sequence(
        ids1=functional_group.get_atom_ids(),
        ids2=clone.get_atom_ids(),
        id_map=id_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_placer_ids(),
        ids2=clone.get_placer_ids(),
        id_map=id_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_core_atom_ids(),
        ids2=clone.get_core_atom_ids(),
        id_map=id_map,
    )


def is_modified_id_sequence(ids1, ids2, id_map):
    for id1, id2 in it.zip_longest(ids1, ids2):
        if id1 in id_map:
            assert id2 == id_map[id1]
        else:
            assert id2 == id1
