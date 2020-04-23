import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def dibromo(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _dibromo(stk.Br(a), stk.C(b), stk.Br(c), stk.C(d))


def _dibromo(bromine1, atom1, bromine2, atom2):
    bonders = (atom1, atom2)
    deleters = (bromine1, bromine2)
    return GenericCaseData(
        functional_group=stk.Dibromo(
            bromine1=bromine1,
            atom1=atom1,
            bromine2=bromine2,
            atom2=atom2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(bromine1, atom1, bromine2, atom2),
        bonders=bonders,
        deleters=deleters,
    )
