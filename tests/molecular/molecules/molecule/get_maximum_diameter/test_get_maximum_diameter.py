import numpy as np

from ...utilities import get_num_atom_ids, normalize_ids
from .utilities import get_maximum_diameter


def test_get_maximum_diameter(case_data, get_atom_ids):
    """
    Test :meth:`.Molecule.get_maximum_diameter`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the correct atomic
        positions of it atoms.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for
        :meth:`.Molecule.get_maximum_diameter`. This allows the
        testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    Notes
    -----
    This test compares the result of
    :meth:`.Molecule.get_maximum_diameter` to the result of
    :func:`.get_maximum_diameter`, which is a utility function defined
    for the purposes of this test. Because
    :func:`.get_maximum_diameter` is tested independently, in
    :mod:`.test_get_maximum_diameter_helper`, if its tests pass, then
    it can be assumed, that :func:`.get_maximum_diameter` gives
    correct results.

    Now, assuming that :func:`.get_maximum_diameter` passed all of its
    tests, this test compares the results of
    :meth:`.Molecule.get_maximum_diameter` to the results of
    :func:`.get_maximum_diameter`. If the results do not match, the
    fault can be placed on :meth:`.Molecule.get_maximum_diameter`,
    because :func:`.get_maximum_diameter` has already been verified to
    be correct by its own tests.

    """

    _test_get_maximum_diameter(
        molecule=case_data.molecule,
        get_atom_ids=get_atom_ids,
        maximum_diameter=get_maximum_diameter(
            position_matrix=case_data.position_matrix,
            atom_ids=tuple(
                normalize_ids(
                    molecule=case_data.molecule,
                    ids=get_atom_ids(case_data.molecule),
                )
            ),
        ),
    )


def _test_get_maximum_diameter(
    molecule,
    get_atom_ids,
    maximum_diameter,
):
    """
    Test :meth:`.Molecule.get_maximum_diameter`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for
        :meth:`.Molecule.get_maximum_diameter`. This allows the
        testing of different values of this parameter.

    maximum_diameter : :class:`float`
        The correct maximum_diameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    num_atom_ids = get_num_atom_ids(molecule, get_atom_ids)
    if num_atom_ids == 1:
        result = molecule.get_maximum_diameter(
            atom_ids=get_atom_ids(molecule),
        )
        assert result == 0
        return

    assert np.allclose(
        a=maximum_diameter,
        b=molecule.get_maximum_diameter(get_atom_ids(molecule)),
        atol=1e-32,
    )
