def molecule_to_json(molecule):
    """
    Return a JSON representation of `molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to serialize.

    Returns
    -------
    :class:`dict`
        A JSON representation of `molecule`.

    """

    return {
        'atoms': tuple(map(atom_to_json, molecule.get_atoms())),
        'bonds': tuple(map(bond_to_json, molecule.get_bonds())),
    }


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

    return {
        'atomic_number': atom.get_atomic_number(),
        'charge': atom.get_charge(),
    }


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

    return {
        'atom1': bond.get_atom1().get_id(),
        'atom2': bond.get_atom2().get_id(),
        'order': bond.get_order(),
        'periodicity': bond.get_periodicity(),
    }
