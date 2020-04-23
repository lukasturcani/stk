import pytest
import stk

from .case_data import CaseData


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
    The periodicity of a bond along the z axis.

    """

    return periodicity_value


@pytest.fixture
def periodicity(x_periodicity, y_periodicity, z_periodicity):
    """
    A possible bond periodicity.

    """

    return x_periodicity, y_periodicity, z_periodicity


@pytest.fixture
def case_data(atom1, atom2, order, periodicity):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        bond=stk.Bond(atom1, atom2, order, periodicity),
        atom1=atom1,
        atom2=atom2,
        order=order,
        periodicity=periodicity,
    )


@pytest.fixture
def bond(case_data):
    """
    A :class:`.Bond` instance.

    """

    return case_data.bond
