def test_get_canonical_atom_ids(case_data):
    """
    Test :meth:`.Molecule.get_canonical_atom_ids`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the expected
        canonical atom ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_canonical_atom_ordering(
        molecule=case_data.molecule,
        canonical_atom_ids=case_data.canonical_atom_ids,
    )


def _test_get_canonical_atom_ordering(molecule, canonical_atom_ids):
    assert molecule.get_canonical_atom_ids() == canonical_atom_ids
