import pytest
import stk


@pytest.fixture(
    params=[
        lambda bond: {},
        lambda bond: {bond.get_atom1().get_id(): stk.Be(33)},
        lambda bond: {bond.get_atom2().get_id(): stk.Bi(122)},
        lambda bond: {
            bond.get_atom1().get_id(): stk.K(4),
            bond.get_atom2().get_id(): stk.S(7),
        },
    ],
)
def get_atom_map(request):
    """
    Return a valid `atom_map` parameter for a bond.

    Parameters
    ----------
    bond : :class:`.Bond`
        The bond, for which the `atom_map` parameter needs to be
        created.

    Returns
    -------
    :class:`dict`
        A valid `atom_map` parameter for :meth:`.Bond.with_atoms`.

    """

    return request.param
