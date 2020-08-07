def test_write(case_data):
    """
    Test writing of molecule to a string.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_write(
        molecule=case_data.molecule,
        writer=case_data.writer,
        string=case_data.string,
    )


def _test_write(molecule, writer, string):
    """
    Test that the written string matches expected string.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        Molecule to test.

    writer : :class:`.MolWriter`
        The writer to test.

    string : :class:`str`
        The expected output string.

    Returns
    -------
    None : :class:`NoneType`

    """

    test_string = writer.to_string(molecule)

    assert test_string == string
