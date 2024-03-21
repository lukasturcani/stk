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

    """
    _test_write(
        molecule=case_data.molecule,
        writer=case_data.writer,
        string=case_data.string,
        periodic_info=case_data.periodic_info,
        file_path=tmp_path / "tmp.coord",
    )


def _test_write(
    molecule: stk.Molecule,
    writer: stk.TurbomoleWriter,
    string: str,
    file_path: pathlib.Path,
    periodic_info: stk.PeriodicInfo | None = None,
):
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

        periodic_info:
            Periodic information.

    """

    writer.write(
        molecule=molecule,
        path=file_path,
        periodic_info=periodic_info,
    )

    with open(file_path, "r") as f:
        content = f.read()

    assert content == string
