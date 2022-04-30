"""
Key Maker Utilities
===================

"""

from __future__ import annotations

import rdkit.Chem.AllChem as rdkit

from ..molecules import Molecule


def get_inchi(
    molecule: Molecule,
):
    """
    Get the InChI of `molecule`.

    Parameters:

        molecule:
            The molecule whose InChI is required.

    Returns:

        The InChI.

    Raises:

        :class:`ValueError`
            If the InChI of `molecule` cannot be generated.

    """

    inchi = rdkit.MolToInchi(molecule.to_rdkit_mol())
    if inchi:
        return inchi
    raise ValueError('The InChI of {molecule} was empty.')


def get_inchi_key(
    molecule: Molecule,
) -> str:
    """
    Get the InChIKey of `molecule`.

    Parameters:

        molecule:
            The molecule whose InChIKey is needed.

    Returns:

        The InChIKey.

    Raises:

        :class:`ValueError`
            If the InChIKey of `molecule` cannot be generated.

    """

    key = rdkit.MolToInchiKey(molecule.to_rdkit_mol())
    if not key:
        raise ValueError(f'InChIKey of {molecule} is empty.')
    return key


def get_smiles(molecule: Molecule) -> str:
    """
    Get the RDKit canonical, isomeric SMILES of `molecule`.

    Parameters:

        molecule:
            The molecule whose SMILES is required.

    Returns:

        The SMILES.

    """

    rdkit_mol = molecule.with_canonical_atom_ordering().to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    rdkit.AssignStereochemistryFrom3D(rdkit_mol)
    rdkit_mol = rdkit.RemoveHs(rdkit_mol)

    return rdkit.MolToSmiles(
        mol=rdkit_mol,
        isomericSmiles=True,
        canonical=True,
    )
