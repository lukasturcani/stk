import numpy as np
import rdkit.Chem.AllChem as rdkit

from ...utilities import is_clone


def test_with_canonical_atom_ordering(case_data):
    """
    Test :meth:`.Molecule.with_canonical_atom_ordering`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the expected
        result.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_with_canonical_atom_ordering(
        molecule=case_data.molecule,
        result=case_data.result,
    )


def _test_with_canonical_atom_ordering(molecule, result):
    """
    Test :meth:`.Molecule.with_canonical_atom_ordering`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    result : :class:`.Molecule`
        The molecule with canonical ordering.

    Returns
    -------
    None : :class:`NoneType`

    """

    ordered = molecule.with_canonical_atom_ordering()
    is_clone(ordered, result)
    order = rdkit.CanonicalRankAtoms(molecule.to_rdkit_mol())
    old_position_matrix = molecule.get_position_matrix()
    new_position_matrix = ordered.get_position_matrix()
    for old_id, new_id in enumerate(order):
        assert np.all(
            np.equal(
                old_position_matrix[old_id],
                new_position_matrix[new_id],
            )
        )
