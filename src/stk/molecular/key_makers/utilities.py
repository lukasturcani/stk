import rdkit.Chem.AllChem as rdkit


def get_inchi(molecule):
    """
    Get the InChI of `molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule whose InChI is required.

    Returns
    -------
    :class:`str`
        The InChI.

    """

    return rdkit.MolToInchi(molecule.to_rdkit_mol())


def get_inchi_key(molecule):
    """
    Get the InChIKey of `molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule whose InChIKey is needed.

    Returns
    -------
    :class:`str`
        The InChIKey.

    """

    return rdkit.MolToInchiKey(molecule.to_rdkit_mol())
