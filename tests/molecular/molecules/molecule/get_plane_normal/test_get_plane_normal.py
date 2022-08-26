import numpy as np

from ...utilities import get_num_atom_ids, normalize_ids
from ..utilities import get_direction
from .utilities import get_plane_normal


def test_get_plane_normal(case_data, get_atom_ids):
    """
    Test :meth:`.Molecule.get_plane_normal`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the correct
        positions of its atoms.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.get_plane_normal`.
        This allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    Notes
    -----
    This test compares the result of
    :meth:`.Molecule.get_plane_normal` to the result of
    :func:`.get_plane_normal`, which is a utility function defined
    for the purposes of this test. Because
    :func:`.get_plane_normal` is tested independently, in
    :mod:`.test_get_plane_normal_helper`, if its tests pass, then
    it can be assumed, that :func:`.get_plane_normal` gives
    correct results.

    Now, assuming that :func:`.get_plane_normal` passed all of its
    tests, this test compares the results of
    :meth:`.Molecule.get_plane_normal` to the results of
    :func:`.get_plane_normal`. If the results do not match, the
    fault can be placed on :meth:`.Molecule.get_plane_normal`,
    because :func:`.get_plane_normal` has already been verified to
    be correct by its own tests.

    """

    _test_get_plane_normal(
        molecule=case_data.molecule,
        get_atom_ids=get_atom_ids,
        normal=get_plane_normal(
            position_matrix=case_data.position_matrix,
            atom_ids=tuple(
                normalize_ids(
                    molecule=case_data.molecule,
                    ids=get_atom_ids(case_data.molecule),
                )
            ),
        ),
    )


def _test_get_plane_normal(molecule, get_atom_ids, normal):
    """
    Test :meth:`.Molecule.get_plane_normal`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.get_plane_normal`.
        This allows the testing of different values of this parameter.

    normal : :class:`numpy.ndarray`
        The correct plane normal.

    Returns
    -------
    None : :class:`NoneType`

    """

    num_atom_ids = get_num_atom_ids(molecule, get_atom_ids)
    if num_atom_ids == 1:
        # Any non-0 vector is valid in this case.
        assert not np.allclose(
            a=[0, 0, 0],
            b=molecule.get_plane_normal(get_atom_ids(molecule)),
            atol=1e-8,
        )
        return

    elif num_atom_ids == 2:
        # Any perpendicular vector is valid in this case.
        result = molecule.get_plane_normal(get_atom_ids(molecule))
        direction = get_direction(
            position_matrix=molecule.get_position_matrix(),
            atom_ids=(0, 1),
        )
        assert result @ direction < 1e-8
        return

    result = molecule.get_plane_normal(get_atom_ids(molecule))
    # The normal may be parallel or anti-parallel.
    assert abs(abs(result @ normal) - 1) < 1e-8
