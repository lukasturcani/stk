import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def thioacid(get_atom_ids):
    a, b, c, d, e = get_atom_ids(5)
    return _thioacid(stk.C(a), stk.O(b), stk.S(c), stk.H(d), stk.C(e))


def _thioacid(carbon, oxygen, sulfur, hydrogen, atom):
    bonders = ()
    deleters = ()
    return GenericCaseData(
        functional_group=stk.Thioacid(
            carbon=carbon,
            oxygen=oxygen,
            sulfur=sulfur,
            hydrogen=hydrogen,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen, sulfur, hydrogen, atom),
        bonders=bonders,
        deleters=deleters,
    )
