from itertools import chain


def _is_metal_atom(atom):
    """
    Check if atom has atomic number of a metal atom.

    Parameters
    ----------
    atom : :class:`.Atom`
        An stk Atom.

    Returns
    -------
    :class:`bool`
        ``True`` if atom has the atomic number of a metal atom.

    """

    # Metal atomic numbers.
    metal_atomic_numbers = set(chain(
        list(range(21, 31)),
        list(range(39, 49)),
        list(range(72, 81))
    ))

    return atom.get_atomic_number() in metal_atomic_numbers
