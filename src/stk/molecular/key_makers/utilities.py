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

    return rdkit.MolToInchi(
        molecule.to_rdkit_mol(),
        treatWarningAsError=True,
    )


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

    """

    rdkit_mol = rdkit.RemoveHs(molecule.to_rdkit_mol())
    return rdkit.MolToSmiles(
        mol=rdkit_mol,
        isomericSmiles=True,
        canonical=True,
    )
