import rdkit.Chem.AllChem as rdkit


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
