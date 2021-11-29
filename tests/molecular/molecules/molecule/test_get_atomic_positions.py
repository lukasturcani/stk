import itertools as it

import numpy as np

from ..utilities import normalize_ids


def test_get_atomic_positions(case_data, get_atom_ids):
    """
    Test :meth:`.Molecule.get_atomic_positions`.

    Parameters
    ----------
    case_data : :class:`.CaseData`.
        A test case. Holds the molecule to test and the correct
        atomic positions.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule` and returns a valid
        `atom_ids` parameter for
        :meth:`.Molecule.get_atomic_positions`. This allows the testing
        of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_atomic_positions(
        molecule=case_data.molecule,
        position_matrix=case_data.position_matrix,
        get_atom_ids=get_atom_ids,
    )


def _test_get_atomic_positions(
    molecule,
    position_matrix,
    get_atom_ids,
):
    """
    Test :meth:`.Molecule.get_atomic_positions`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    position_matrix : :class:`numpy.ndarray`
        The correct atomic positions.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule` and returns a valid
        `atom_ids` parameter for
        :meth:`.Molecule.get_atomic_positions`. This allows the testing
        of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    for atom_id, position in it.zip_longest(
        normalize_ids(molecule, get_atom_ids(molecule)),
        molecule.get_atomic_positions(get_atom_ids(molecule)),
    ):
        assert np.allclose(
            a=position,
            b=position_matrix[atom_id],
            atol=1e-6,
        )
