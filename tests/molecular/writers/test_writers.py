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
        periodic_cell=case_data.periodic_cell,
    )


def _test_write(molecule, writer, string, periodic_cell=None):
    """
    Test that the written string matches expected string.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        Molecule to test.

    writer : :class:`.???????`
        The writer to test.

    string : :class:`str`
        The expected output string.

    periodic_cell : :class:`tuple` of :class:`np.array`, optional
        Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

    Returns
    -------
    None : :class:`NoneType`

    """

    test_string = writer.write(
        mol=molecule,
        periodic_cell=periodic_cell,
    )

    assert test_string == string
