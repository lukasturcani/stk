import stk

from .case_data import CaseData


def test_to_string(case_data: CaseData) -> None:
    """
    Test writing of molecule to a string.

    Parameters:
        case_data:
            A test case.

    """

    _test_to_string(
        molecule=case_data.molecule,
        writer=case_data.writer,
        string=case_data.string,
    )


def _test_to_string(
    molecule: stk.Molecule,
    writer: stk.XyzWriter,
    string: str,
) -> None:
    """
    Test that the written string matches expected string.

    Parameters:
        molecule:
            Molecule to test.

        writer:
            The writer to test.

        string:
            The expected output string.

    """

    test_string = writer.to_string(molecule)

    assert test_string == string
