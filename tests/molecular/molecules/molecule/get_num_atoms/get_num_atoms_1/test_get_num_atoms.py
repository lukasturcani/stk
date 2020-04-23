def test_get_num_atoms(case_data):
    """
    Test :meth:`.Molecule.get_num_atoms`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the number of atoms
        it should have.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.molecule.get_num_atoms() == case_data.num_atoms
