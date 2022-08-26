import numpy as np

from ....utilities import is_clone


def test_with_structure_from_file(datadir, case_data):
    """
    Test :meth:`.Molecule.with_structure_from_file`.

    Parameters
    ----------
    datadir : :class:`pathlib.Path`
        Holds the data files for the test.

    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the path of its
        structure file.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save a copy of the position matrix, to ensure the original
    # molecule is not modified by the test, because it is meant to be
    # immutable.
    position_matrix = case_data.molecule.get_position_matrix()
    _test_with_structure_from_file(
        molecule=case_data.molecule,
        path=str(datadir / case_data.path),
    )
    assert np.all(
        np.equal(
            position_matrix,
            case_data.molecule.get_position_matrix(),
        )
    )


def _test_with_structure_from_file(molecule, path):
    """
    Test :meth:`.Molecule.with_structure_from_file`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    path : :class:`str`
        The path the structure file of `molecule`.

    Returns
    -------
    None : :class:`NoneType`

    """

    new = molecule.with_structure_from_file(path)
    is_clone(new, molecule)
    size_diff = abs(
        molecule.get_maximum_diameter() - new.get_maximum_diameter()
    )
    assert size_diff > 1
