import numpy as np

from tests.utilities import is_equivalent_molecule


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
        molecule=case_data.molecule,
    )


def _test_from_json(
    dejsonizer,
    json,
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

    molecule : :class:`.Molecule`
        The correct dejsonized molecule.

    Returns
    -------
    None : :class:`NoneType`

    """

    result = dejsonizer.from_json(json)
    is_equivalent_molecule(result, molecule)
    assert np.all(np.equal(
        molecule.get_position_matrix(),
        result.get_position_matrix(),
    ))
