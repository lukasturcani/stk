import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def diol(get_atom_ids):
    a, b, c, d, e, f = get_atom_ids(6)
    return _diol(
        atom1=stk.C(a),
        oxygen1=stk.O(b),
        hydrogen1=stk.H(c),
        atom2=stk.C(d),
        oxygen2=stk.O(e),
        hydrogen2=stk.H(f),
    )


def _diol(atom1, oxygen1, hydrogen1, atom2, oxygen2, hydrogen2):
    bonders = (atom1, atom2)
    deleters = (oxygen1, hydrogen1, oxygen2, hydrogen2)
    return GenericCaseData(
        functional_group=stk.Diol(
            atom1=atom1,
            oxygen1=oxygen1,
            hydrogen1=hydrogen1,
            atom2=atom2,
            oxygen2=oxygen2,
            hydrogen2=hydrogen2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(atom1, oxygen1, hydrogen1, atom2, oxygen2, hydrogen2),
        bonders=bonders,
        deleters=deleters,
    )
