from typing import Callable

import pytest
import stk


@pytest.fixture(
    params=[
        lambda bond: {},
        lambda bond: {bond.get_atom1().get_id(): 33},
        lambda bond: {bond.get_atom2().get_id(): 122},
        lambda bond: {
            bond.get_atom1().get_id(): 4,
            bond.get_atom2().get_id(): 7,
        },
    ],
)
def get_id_map(request) -> Callable[[stk.Bond], dict[int, int]]:
    """
    Return a valid `atom_map` parameter for a bond.

    Parameters:

        bond:
            The bond, for which the `atom_map` parameter needs to be
            created.

    Returns:

        A valid `id_map` parameter for :meth:`.Bond.id_atoms`.

    """

    return request.param
