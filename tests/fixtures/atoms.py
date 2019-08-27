import pytest
import stk


@pytest.fixture('session')
def hydrogen():
    return stk.H(1, -3, attr1='12', attr2=32, _attr3=12.2)


@pytest.fixture('session')
def carbon():
    return stk.C(333, attr2=None, attr12='abc', _pattr=22.22)


@pytest.fixture('session')
def lithium():
    return stk.Li(121, 2, alpha=232, beta='bvc', _gamma=4523)


@pytest.fixture('session')
def chlorine():
    return stk.Cl(786, a=9, b='bbb', _c='yfg')
