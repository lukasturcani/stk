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
    A function which takes a tuple and returns a subset of it.

    """

    return request.param


@pytest.fixture
def get_bonders(get_subset):
    """
    A function which takes atoms and returns the bonders.

    """

    return get_subset


@pytest.fixture
def get_deleters(get_subset):
    """
    A function which takes atoms and returns the deleters.

    """

    return get_subset


@pytest.fixture(
    params=[
        stk.Amine,
        stk.Aldehyde,
        stk.CarboxylicAcid,
        stk.Amide,
        stk.Thioacid,
        stk.Alcohol,
        stk.Thiol,
        stk.Fluoro,
        stk.Bromo,
        stk.Iodo,
        stk.TerminalAlkyne,
        stk.TerminalAlkene,
        stk.BoronicAcid,
        stk.Diol,
        stk.Difluoro,
        stk.Dibromo,
        stk.RingAmine,
    ],
)
def get_functional_group(request):
    """
    A function which creates a :class:`.FunctionalGroup` instance.

    The function must take atoms, bonders and deleters parameters,
    which are the atoms the :class:`.FunctionalGroup` should hold.

    """

    return request.param
