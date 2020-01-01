import pytest
import stk


@pytest.fixture(
    params=[
            (stk.N(0), stk.H(21), stk.K(3), stk.Li(223), ),
            (stk.H(120), ),
            (),
    ],
)
def atoms(request):
    """
    A :class:`tuple` of :class:`.Atom` instances.

    """

    return tuple(a.clone() for a in request.param)


@pytest.fixture(
    params=[
        lambda atoms: tuple(atoms),
        lambda atoms: (),
        lambda atoms: atoms[0:1] if atoms else (),
        lambda atoms: tuple(atoms[i] for i in range(0, len(atoms), 2))
    ],
)
def get_subset(request):
    """
    Return a subset of a :class:`tuple`.

    Parameters
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        A collection of atoms, from which a subset is returned.

    Returns
    -------
    :class:`tuple` of :class:`.Atom`
        A subset of `atoms`.

    """

    return request.param


@pytest.fixture
def get_bonders(get_subset):
    """
    Return the subset of `atoms` which are bonders.

    Parameters
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        A collection of atoms, from which the bonders should be
        selected.

    Returns
    -------
    :class:`tuple` of :class:`.Atom`
        The bonder atoms.

    """

    return get_subset


@pytest.fixture
def get_deleters(get_subset):
    """
    Return the subset of `atoms` which are deleters.

    Parameters
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        A collection of atoms, from which the deleters should be
        selected.

    Returns
    -------
    :class:`tuple` of :class:`.Atom`
        The deleter atoms.


    """

    return get_subset


@pytest.fixture(
    params=[
        stk.PrimaryAmine,
        stk.SecondaryAmine,
        stk.Aldehyde,
        stk.CarboxylicAcid,
        stk.Amide,
        stk.Thioacid,
        stk.Alcohol,
        stk.Thiol,
        stk.Fluoro,
        stk.Bromo,
        stk.Iodo,
        stk.Alkyne,
        stk.Alkene,
        stk.BoronicAcid,
        stk.Diol,
        stk.Difluoro,
        stk.Dibromo,
        stk.RingAmine,
    ],
)
def get_functional_group(request):
    """
    Return a :class:`.FunctionalGroup` instance.

    Parameters
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        The atoms in the functional group.

    bonders : :class:`tuple` of :class:.Atom`
        The bonder atoms in the functional group.

    deleters : :class:`tuple` of :class:.Atom`
        The deleter atoms in the functional group.

    Returns
    -------
    :class:`.FunctionalGroup`
        A functional group.

    """

    return request.param
