"""
MAE Updating Utilities
======================

"""


from stk.utilities import mol_from_mae_file


def _with_structure_from_mae(self, path):
    """
    Change structure to match an ``.mae`` file.

    Parameters
    ----------
    path : :class:`str`
        The full path of the ``.mae`` file from which the structure
        should be updated.

    Returns
    -------
    :class:`.Molecule`
        The molecule.

    """

    molecule = mol_from_mae_file(path)
    return self._with_position_matrix(
        position_matrix=molecule.GetConformer().GetPositions()
    )
