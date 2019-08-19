import pytest
import stk


@pytest.fixture('session')
def bond(hydrogen, carbon):
    return stk.Bond(
        atom1=hydrogen,
        atom2=carbon,
        order=2,
        attr1=1,
        attr2='2',
        _attr3=12.2
    )


@pytest.fixture('session')
def periodic_bond(lithium, chlorine):
    return stk.Bond(
        atom1=lithium,
        atom2=chlorine,
        order=21,
        periodicity=(1, 0, -1),
        attr10=16,
        attr20='26',
        _attr30=126.2
    )
