import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def aldehyde(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _aldehyde(stk.C(a), stk.O(b), stk.H(c), stk.C(d))


def _aldehyde(carbon, oxygen, hydrogen, atom):
    bonders = (carbon,)
    deleters = (oxygen,)
    return GenericCaseData(
        functional_group=stk.Aldehyde(
            carbon=carbon,
            oxygen=oxygen,
            hydrogen=hydrogen,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen, hydrogen, atom),
        bonders=bonders,
        deleters=deleters,
    )
