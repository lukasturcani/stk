import numpy as np


def test_get_position_matrix(case_data):
    """
    Test :meth:`.Molecule.get_position_matrix`.

    Parameters
    ----------
    case_data : :class:`.CaseData
        A test case. Holds the molecule to test and the correct
        position matrix.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.allclose(
        a=case_data.position_matrix,
        b=case_data.molecule.get_position_matrix(),
        atol=1e-6,
    )
