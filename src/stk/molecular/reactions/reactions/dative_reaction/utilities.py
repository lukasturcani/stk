"""
Dative Reaction Utilities
=========================

"""

from itertools import chain


def is_metal(atom):
    """
    Check if `atom` is a metal atom.

    Parameters
    ----------
    atom : :class:`.Atom`
        An atom.

    Returns
    -------
    :class:`bool`
        ``True`` if `atom` is a metal atom.

    """

    metal_atomic_numbers = chain(
        range(21, 31),
        range(39, 49),
        range(72, 81),
    )
    return atom.get_atomic_number() in metal_atomic_numbers
