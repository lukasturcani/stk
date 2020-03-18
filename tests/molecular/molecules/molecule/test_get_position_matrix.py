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

    assert np.all(np.equal(
        case_data.position_matrix,
        case_data.molecule.get_position_matrix(),
    ))
