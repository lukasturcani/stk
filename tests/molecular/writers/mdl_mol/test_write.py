import stk

from .case_data import CaseData
import pathlib


def test_write(case_data: CaseData, tmp_path: pathlib.Path) -> None:
    """
    Test writing of molecule to a file.

    Parameters:
        case_data:
            A test case.

        tmp_path:
            Path to temporary directory.

    Returns:
        :class:`NoneType`

    """

    _test_write(
        molecule=case_data.molecule,
        writer=case_data.writer,
        string=case_data.string,
        file_path=tmp_path / "tmp.mol",
    )


def _test_write(
    molecule: stk.Molecule,
    writer: stk.MolWriter,
    string: str,
    file_path: pathlib.Path,
) -> None:
    """
    Test that the written file content matches expected string.

    Parameters:
        molecule:
            Molecule to test.

        writer:
            The writer to test.

        string:
            The expected output string.

        file_path:
            Path to temporary file.

    Returns:
        :class:`NoneType`

    """

    writer.write(molecule, file_path)

    with open(file_path, "r") as f:
        content = f.read()

    assert content == string
