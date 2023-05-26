import pytest
import stk


def get_atom_map_0(functional_group):
    """
    Get an atom_map with a single atom.

    """

    # Make sure new_id is always valid, by making 1 larger than the
    # biggest on in the functional group. This prevents two atoms in
    # the functional group from having the same id.
    new_id = max(functional_group.get_atom_ids()) + 1
    atoms = (stk.Li(new_id),)
    return dict(zip(functional_group.get_atom_ids(), atoms))


@pytest.fixture(
    params=[
        lambda functional_group: {},
        get_atom_map_0,
    ],
)
def get_atom_map(request):
    """
    Return a valid `atom_map` parameter for a functional group.

    Parameters
    ----------
    functional_group : :class:`.FunctionalGroup`
        A functional group being cloned, for which a valid
        `atom_map` parameter needs to created.

    Returns
    -------
    :class:`dict`
        A valid `atom_map` parameter for calling
        :meth:`~.FunctionalGroup.with_atoms` on `functional_group`.

    """

    return request.param
