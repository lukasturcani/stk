import itertools as it
from typing import Callable, Iterable

import stk

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
        original_ids=functional_group.get_atom_ids(),
        modified_ids=clone.get_atom_ids(),
        id_map=id_map,
    )
    is_modified_id_sequence(
        original_ids=functional_group.get_placer_ids(),
        modified_ids=clone.get_placer_ids(),
        id_map=id_map,
    )
    is_modified_id_sequence(
        original_ids=functional_group.get_core_atom_ids(),
        modified_ids=clone.get_core_atom_ids(),
        id_map=id_map,
    )


def is_modified_id_sequence(
    original_ids: Iterable[int],
    modified_ids: Iterable[int],
    id_map: dict[int, int],
) -> None:
    for original_id, modified_id in it.zip_longest(
        original_ids,
        modified_ids,
    ):
        if original_id in id_map:
            assert modified_id == id_map[original_id]
        else:
            assert original_id == modified_id
