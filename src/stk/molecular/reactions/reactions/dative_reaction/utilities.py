"""
Dative Reaction Utilities
=========================

"""

from itertools import chain

from ....atom import Atom

__all__ = (
    'is_metal',
)


def is_metal(atom: Atom) -> bool:
    """
    Check if `atom` is a metal atom.

    Parameters:

        atom:
            An atom.

    Returns:

        ``True`` if `atom` is a metal atom.

    """

    metal_atomic_numbers = chain(
        range(21, 31),
        range(39, 49),
        range(72, 81),
    )
    return atom.get_atomic_number() in metal_atomic_numbers
