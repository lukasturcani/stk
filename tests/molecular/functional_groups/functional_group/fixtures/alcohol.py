import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def alcohol(get_atom_ids):
    a, b, c = get_atom_ids(3)
    return _alcohol(stk.O(a), stk.H(b), stk.C(c))


def _alcohol(oxygen, hydrogen, atom):
    bonders = (oxygen,)
    deleters = (hydrogen,)
    return GenericCaseData(
        functional_group=stk.Alcohol(
            oxygen=oxygen,
            hydrogen=hydrogen,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(oxygen, hydrogen, atom),
        bonders=bonders,
        deleters=deleters,
    )
