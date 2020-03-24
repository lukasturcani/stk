import numpy as np

from tests.utilities import is_equivalent_constructed_molecule


def test_from_json(case_data):
    """
    Test :meth:`.ConstructedMoleculeDejsonizer.from_json`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the dejsonizer to test and the correct
        dejsonized molecule.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_from_json(
        dejsonizer=case_data.dejsonizer,
        molecule_json=case_data.molecule_json,
        constructed_molecule_json=case_data.constructed_molecule_json,
        position_matrix=case_data.position_matrix,
        building_blocks=case_data.building_blocks,
        molecule=case_data.molecule,
    )


def _test_from_json(
    dejsonizer,
    molecule_json,
    constructed_molecule_json,
    position_matrix,
    building_blocks,
    molecule,
):
    """
    Test :meth:`.ConstructedMoleculeDejsonizer.from_json`.

    Parameters
    ----------
    molecule_json : :class:`dict`
        A JSON of the molecular information of the constructed
        molecule.

    constructed_molecule_json : :class:`dict`
        A JSON of the constructed molecule information of the
        constructed molecule.

    position_matrix : :class:`numpy.ndarray`
        The position matrix of the constructed molecule.

    building_blocks : :class:`tuple` of :class:`.Molecule`
        The building blocks of the constructed molecule.

    molecule : :class:`.ConstructedMolecule`
        The correct dejsonized molecule.

    Returns
    -------
    None : :class:`NoneType`

    """

    result = dejsonizer.from_json(
        molecule_json=molecule_json,
        constructed_molecule_json=constructed_molecule_json,
        position_matrix=position_matrix,
        building_blocks=building_blocks,
    )
    is_equivalent_constructed_molecule(result, molecule)
    assert np.all(np.equal(
        result.get_position_matrix(),
        molecule.get_position_matrix(),
    ))
