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

    Raises
    ------
    :class:`ValueError`
        If InChI of `molecule` is an empty string.

    """

    key = rdkit.MolToInchi(
        mol=molecule.to_rdkit_mol(),
        treatWarningAsError=True,
    )
    if not key:
        raise ValueError(
            f'InChI of {molecule} is empty string'
        )

    return key


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

    Raises
    ------
    :class:`ValueError`
        If InChIKey of `molecule` is an empty string.

    """

    key = rdkit.MolToInchiKey(molecule.to_rdkit_mol())
    if not key:
        raise ValueError(
            f'InChIKey of {molecule} is empty string'
        )

    return key


def get_smiles(molecule):
    """
    Get the RDKit canonical, isomeric SMILES of `molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule whose SMILES is required.

    Returns
    -------
    :class:`str`
        The SMILES.

    Raises
    ------
    :class:`ValueError`
        If SMILES of `molecule` is an empty string.

    """

    rdkit_mol = molecule.to_rdkit_mol()
    rdkit.AssignStereochemistryFrom3D(rdkit_mol)
    rdkit.SanitizeMol(rdkit_mol)
    rdkit_mol = rdkit.RemoveHs(rdkit_mol)
    key = rdkit.MolToSmiles(
        mol=rdkit_mol,
        isomericSmiles=True,
        canonical=True,
    )
    if not key:
        raise ValueError(
            f'SMILES of {molecule} is empty string'
        )

    return key
