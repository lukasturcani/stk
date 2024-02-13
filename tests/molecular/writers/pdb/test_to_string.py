import stk

from .case_data import CaseData


def test_to_string(case_data: CaseData) -> None:
    """
    Test writing of molecule to a string.

    Parameters:
        case_data:
            A test case.

    Returns:
        :class:`NoneType`

    """

    _test_to_string(
        molecule=case_data.molecule,
        writer=case_data.writer,
        string=case_data.string,
        periodic_info=case_data.periodic_info,
    )


def _test_to_string(
    molecule: stk.Molecule,
    writer: stk.PdbWriter,
    string: str,
    periodic_info: stk.PeriodicInfo | None = None,
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

        periodic_info:
            Periodic information.

    Returns:
        :class:`NoneType`

    """

    test_string = writer.to_string(
        molecule=molecule,
        periodic_info=periodic_info,
    )

    assert test_string == string
