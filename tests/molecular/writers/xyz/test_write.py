def test_write(case_data, tmp_path):
    """
    Test writing of molecule to a file.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case.

    tmp_path : :class:`pathlib2.Path`
        Path to temporary directory.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_write(
        molecule=case_data.molecule,
        writer=case_data.writer,
        string=case_data.string,
        file_path=tmp_path / "tmp.xyz",
    )


def _test_write(molecule, writer, string, file_path):
    """
    Test that the written file content matches expected string.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        Molecule to test.

    writer : :class:`.XyzWriter`
        The writer to test.

    string : :class:`str`
        The expected output string.

    file_path : :class:`str`
        Path to temporary file.

    Returns
    -------
    None : :class:`NoneType`

    """

    writer.write(molecule, file_path)

    with open(file_path, "r") as f:
        content = f.read()

    assert content == string
