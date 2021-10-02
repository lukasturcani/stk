"""
Functional Group Factory Utilities
==================================

"""

from collections import abc
import rdkit.Chem.AllChem as rdkit

from ...molecule import Molecule


def get_atom_ids(
    query: str,
    molecule: Molecule,
) -> abc.Iterable[int]:
    """
    Yield the ids of atoms in `molecule` which match `query`.

    Multiple substructures in `molecule` can match `query` and
    therefore each set is yielded as a group.

    Parameters:

        query:
            A SMARTS string used to query atoms.

        molecule:
            A molecule whose atoms should be queried.

    Yields:

        The ids of atoms in `molecule` which match `query`.

    """

    rdkit_mol = molecule.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    yield from rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(query),
    )
