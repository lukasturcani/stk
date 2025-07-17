import typing

import rdkit.Chem.AllChem as rdkit

from stk._internal.utilities.utilities import remake


def _with_structure_from_pdb(self, path: str) -> typing.Self:  # type: ignore[misc]
    """
    Change structure to match a ``.pdb`` file.

    Parameters:

        path:
            The full path of the ``.pdb`` file from which the structure
            should be updated.

    Returns:
        The molecule.

    """

    rdk_mol = rdkit.MolFromPDBFile(
        pdbFileName=path,
        sanitize=False,
        removeHs=False,
    )
    if rdk_mol is None:
        rdkit.MolFromPDBFile(
            molFileName=path,
            sanitize=False,
            removeHs=False,
        )

    molecule = remake(rdk_mol)
    return self._with_position_matrix(
        position_matrix=molecule.GetConformer().GetPositions()
    )
