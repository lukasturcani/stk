import numpy as np

from ..utilities import is_equivalent_molecule


def test_from_json(case_data):
    """
    Test :meth:`.MoleculeDejsonizer.from_json`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the dejsonizer to test and the correct
        dejsonized molecule.

    Results
    -------
    None : :class:`NoneType`

    """

    _test_from_json(
        dejsonizer=case_data.dejsonizer,
        json=case_data.json,
        position_matrix=case_data.position_matrix,
        molecule=case_data.molecule,
    )


def _test_from_json(
    dejsonizer,
    json,
    position_matrix,
    molecule,
):
    """
    Test :meth:`.MoleculeDejsonizer.from_json`.

    Parameters
    ----------
    dejsonizer : :class:`.MoleculeDejsonizer`
        The dejsonizer to test.

    json : :class:`dict`
        The JSON to test.

    position_matrix : :class:`list`
        The position matrix of the dejsonized molecule.

    molecule : :class:`.Molecule`
        The correct dejsonized molecule.

    Returns
    -------
    None : :class:`NoneType`

    """

    result = dejsonizer.from_json(json, position_matrix)
    is_equivalent_molecule(result, molecule)
    assert np.all(np.equal(
        molecule.get_position_matrix(),
        result.get_position_matrix(),
    ))
