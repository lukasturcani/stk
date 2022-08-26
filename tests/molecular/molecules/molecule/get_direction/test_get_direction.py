import numpy as np

from ...utilities import get_num_atom_ids, normalize_ids
from ..utilities import get_direction


def test_get_direction(case_data, get_atom_ids):
    """
    Test :meth:`.Molecule.get_direction`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and its correct
        atomic positions.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.get_direction`. This
        allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    Notes
    -----
    This test compares the result of :meth:`.Molecule.get_direction`
    to the result of :func:`.get_direction`, which is a utility
    function defined for the purposes of this test. Because
    :func:`.get_direction` is tested independently, in
    :mod:`.test_get_direction_helper`, if its tests pass, then
    it can be assumed, that :func:`.get_direction` gives correct
    results.

    Now, assuming that :func:`.get_direction` passed all of its tests,
    this test compares the results of :meth:`.Molecule.get_direction`
    to the results of :func:`.get_direction`. If the results do not
    match, the fault can be placed on :meth:`.Molecule.get_direction`,
    because :func:`.get_direction` has already been verified to be
    correct by its own tests.

    """

    _test_get_direction(
        molecule=case_data.molecule,
        direction=get_direction(
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


def _test_get_direction(molecule, direction, get_atom_ids):
    """
    Test :meth:`.Molecule.get_direction`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    direction : :class:`.Molecule`
        The correct direction of `molecule`.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.get_direction`. This
        allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    if get_num_atom_ids(molecule, get_atom_ids) == 1:
        # Any non-0 vector is valid in this case.
        assert not np.allclose(
            a=[0, 0, 0],
            b=molecule.get_direction(get_atom_ids(molecule)),
            atol=1e-13,
        )
        return

    result = molecule.get_direction(get_atom_ids(molecule))
    # The direction may be parallel or anti-parallel.
    return abs(abs(result @ direction) - 1) < 1e-13
