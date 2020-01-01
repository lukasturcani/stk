import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture
def primary_amine(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _primary_amine(stk.N(a), stk.H(b), stk.H(c), stk.C(d))


def _primary_amine(nitrogen, hydrogen1, hydrogen2, atom):
    bonders = (nitrogen, )
    deleters = (hydrogen1, hydrogen2)
    return _TestCase(
        functional_group=stk.PrimaryAmine(
            nitrogen=nitrogen,
            hydrogen1=hydrogen1,
            hydrogen2=hydrogen2,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(nitrogen, hydrogen1, hydrogen2, atom),
        bonders=bonders,
        deleters=deleters,
    )
