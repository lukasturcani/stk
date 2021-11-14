from __future__ import annotations

from typing import Callable

import stk


def test_with_ids(
    bond: stk.Bond,
    get_id_map: Callable[[stk.Bond], dict[int, int]],
) -> None:
    """
    Test :meth:`.Bond.with_ids`.

    Parameters:

        bond : :class:`.Bond`
            The bond to test.

        get_id_map:
            Takes a single parameter, `bond`, and returns a valid
            `id_map` parameter for its :meth:`.Bond.with_ids`
            method. This allows the testing of different values of
            this parameter.

    """

    id_map = get_id_map(bond)
    clone = bond.with_ids(id_map)
    assert clone is not bond

    expected_id1 = id_map.get(
        bond.get_atom1().get_id(),
        bond.get_atom1().get_id(),
    )
    assert clone.get_atom1().get_id() == expected_id1

    expected_id2 = id_map.get(
        bond.get_atom2().get_id(),
        bond.get_atom2().get_id(),
    )
    assert clone.get_atom2().get_id() == expected_id2

    assert bond.get_periodicity() == clone.get_periodicity()
    assert bond.get_order() == clone.get_order()
