import pytest
import stk


@pytest.fixture(params=[stk.C(0), stk.H(1), stk.K(12)])
def atom(request):
    """
    An :class:`.Atom` instance.

    """

    return request.param.clone()


@pytest.fixture
def atom1(atom):
    """
    An :class:`.Atom` instance.

    """

    return atom.clone()


@pytest.fixture
def atom2(atom):
    """
    An :class:`.Atom` instance.

    """

    return atom.clone()


@pytest.fixture(params=[1, 2, 3, 4, 1.5])
def order(request):
    """
    A possible bond order.

    """

    return request.param


@pytest.fixture(params=list(range(-3, 4)))
def periodicity_value(request):
    """
    The periodicity of a bond along 1 axis.

    """

    return request.param


@pytest.fixture
def x_periodicity(periodicity_value):
    """
    The periodicity of a bond along the x axis.

    """

    return periodicity_value


@pytest.fixture
def y_periodicity(periodicity_value):
    """
    The periodicity of a bond along the y axis.

    """

    return periodicity_value


@pytest.fixture
def z_periodicity(periodicity_value):
    """
    THe periodicity of a bond along the z axis.

    """

    return periodicity_value


@pytest.fixture
def periodicity(x_periodicity, y_periodicity, z_periodicity):
    """
    A possible bond periodicity.

    """

    return x_periodicity, y_periodicity, z_periodicity


@pytest.fixture
def bond(atom1, atom2, order, periodicity):
    """
    A :class:`.Bond` instance.

    """

    return stk.Bond(atom1, atom2, order, periodicity)
