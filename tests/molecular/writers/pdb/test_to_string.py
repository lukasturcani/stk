def test_to_string(case_data):
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

    _test_to_string(
        molecule=case_data.molecule,
        writer=case_data.writer,
        string=case_data.string,
        periodic_info=case_data.periodic_info,
    )


def _test_to_string(molecule, writer, string, periodic_info=None):
    """
    Test that the written string matches expected string.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        Molecule to test.

    writer : :class:`.PdbWriter`
        The writer to test.

    string : :class:`str`
        The expected output string.

    periodic_info : :class:`.PeriodicInfo`
        Periodic information.

    Returns
    -------
    None : :class:`NoneType`

    """

    test_string = writer.to_string(
        molecule=molecule,
        periodic_info=periodic_info,
    )

    assert test_string == string
