from typing import Callable

import pytest
import stk


def get_id_map_0(
    functional_group: stk.FunctionalGroup,
) -> dict[int, int]:
    """
    Get an id_map with a single atom.

    """

    # Make sure new_id is always valid, by making 1 larger than the
    # biggest on in the functional group. This prevents two atoms in
    # the functional group from having the same id.
    new_id = max(functional_group.get_atom_ids()) + 1
    ids = (new_id,)
    return dict(zip(functional_group.get_atom_ids(), ids))


def get_id_map_1(
    functional_group: stk.FunctionalGroup,
) -> dict[int, int]:
    new_id = max(functional_group.get_atom_ids()) + 1
    num_atoms = sum(1 for _ in functional_group.get_atoms())
    ids = range(new_id, new_id + num_atoms)
    return dict(zip(functional_group.get_atom_ids(), ids))


@pytest.fixture(
    params=[
        lambda functional_group: {},
        get_id_map_0,
        get_id_map_1,
    ],
)
def get_id_map(
    request,
) -> Callable[[stk.FunctionalGroup], dict[int, int]]:
    """
    Return a valid `id_map` parameter for a functional group.

    Parameters:

        functional_group:
            A functional group being tested, for which a valid
            `id_map` parameter needs to created.

    Returns:

        A valid `id_map` parameter for calling
        :meth:`~.FunctionalGroup.with_ids` on `functional_group`.

    """

    return request.param
