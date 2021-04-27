"""
Functional Group Factory Utilities
==================================

"""

import rdkit.Chem.AllChem as rdkit


def _get_atom_ids(query, molecule):
    """
    Yield the ids of atoms in `molecule` which match `query`.

    Multiple substructures in `molecule` can match `query` and
    therefore each set is yielded as a group.

    Parameters
    ----------
    query : :class:`str`
        A SMARTS string used to query atoms.

    molecule : :class:`.Molecule`
        A molecule whose atoms should be queried.

    Yields
    ------
    :class:`tuple` of :class:`int`
        The ids of atoms in `molecule` which match `query`.

    """

    rdkit_mol = molecule.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    yield from rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(query),
    )
