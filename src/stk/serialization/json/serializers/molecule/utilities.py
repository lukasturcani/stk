"""
JSON Serialization Utilities for Molecules
==========================================

"""


def atom_to_json(atom):
    """
    Return a JSON representation of `atom`.

    Parameters
    ----------
    atom : :class:`.Atom`
        The atom to serialize.

    Returns
    -------
    :class:`dict`
        A JSON representation of `atom`.

    """

    return (
        atom.get_atomic_number(),
        atom.get_charge(),
    )


def bond_to_json(bond):
    """
    Return a JSON representation of `bond`.

    Parameters
    ----------
    bond : :class:`.Bond`
        The bond to serialize.

    Returns
    -------
    :class:`dict`
        A JSON representation of `bond`.

    """

    return (
        bond.get_atom1().get_id(),
        bond.get_atom2().get_id(),
        bond.get_order(),
        bond.get_periodicity(),
    )
