import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture
def secondary_amine(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _secondary_amine(stk.N(a), stk.H(b), stk.C(c), stk.C(d))


def _secondary_amine(nitrogen, hydrogen, atom1, atom2):
    bonders = (nitrogen, )
    deleters = (hydrogen, )
    return _TestCase(
        functional_group=stk.SecondaryAmine(
            nitrogen=nitrogen,
            hydrogen=hydrogen,
            atom1=atom1,
            atom2=atom2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(nitrogen, hydrogen, atom1, atom2),
        bonders=bonders,
        deleters=deleters,
    )
