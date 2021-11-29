"""
MongoDB Utilities
=================

"""

from collections.abc import Sequence
from typing import Iterable, Optional


class HashableDict(dict):
    def __hash__(self):
        return hash((
            frozenset(self),
            frozenset(_to_hashable(list(self.values()))),
        ))

    def __eq__(self, other):
        return super().__eq__(other)


def _to_hashable(item):
    if isinstance(item, list):
        return tuple(_to_hashable(subitem) for subitem in item)
    else:
        return item


def get_any_value(
    mapping: dict[str, Sequence[dict]],
    keys: Iterable[str],
) -> Optional[dict]:
    """
    Return any value in `mapping` for any of the `keys`.

    This function will only return a value if it is a non-empty
    sequence. The function will not return the sequence itself, but
    its first member.

    Parameters:

        mapping:
            A mapping from which a value is to be extracted.

        keys:
            The keys which are used to look up values in `mapping`.

    Returns:

        The value of the first key found in `mapping` wll be returned.
        If none of the `keys` are found in `mapping`, ``None`` will be
        returned.

    """

    for key in keys:
        if key in mapping and len(mapping[key]) > 0:
            return mapping[key][0]
    return None
