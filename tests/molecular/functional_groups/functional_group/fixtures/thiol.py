import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture
def thiol(get_atom_ids):
    a, b, c = get_atom_ids(3)
    return _thiol(stk.S(a), stk.H(b), stk.C(c))


def _thiol(sulfur, hydrogen, atom):
    bonders = ()
    deleters = ()
    return _TestCase(
        functional_group=stk.Thiol(
            sulfur=sulfur,
            hydrogen=hydrogen,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(sulfur, hydrogen, atom),
        bonders=bonders,
        deleters=deleters,
    )
