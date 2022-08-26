import numpy as np

from ...utilities import normalize_ids
from ..utilities import get_centroid


def test_get_centroid(case_data, get_atom_ids):
    """
    Test :meth:`.Molecule.get_centroid`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the correct atom
        positions.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.get_centroid`. This
        allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    Notes
    -----
    This test compares the result of :meth:`.Molecule.get_centroid`
    to the result of :func:`.get_centroid`, which is a utility
    function defined for the purposes of this test. Because
    :func:`.get_centroid` is tested independently, in
    :mod:`.test_get_centroid_helper`, if its tests pass, then
    it can be assumed, that :func:`.get_centroid` gives correct
    results.

    Now, assuming that :func:`.get_centroid` passed all of its tests,
    this test compares the results of :meth:`.Molecule.get_centroid`
    to the results of :func:`.get_centroid`. If the results do not
    match, the fault can be placed on :meth:`.Molecule.get_centroid`,
    because :func:`.get_centroid` has already been verified to be
    correct by its own tests.

    """

    _test_get_centroid(
        molecule=case_data.molecule,
        centroid=get_centroid(
            position_matrix=case_data.position_matrix,
            atom_ids=tuple(
                normalize_ids(
                    molecule=case_data.molecule,
                    ids=get_atom_ids(case_data.molecule),
                )
            ),
        ),
        get_atom_ids=get_atom_ids,
    )


def _test_get_centroid(molecule, centroid, get_atom_ids):
    """
    Test :meth:`.Molecule.get_centroid`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    centroid : :class:`numpy.ndarray`
        The correct centroid of `molecule`.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.get_centroid`. This
        allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.allclose(
        a=centroid,
        b=molecule.get_centroid(get_atom_ids(molecule)),
        atol=1e-6,
    )
