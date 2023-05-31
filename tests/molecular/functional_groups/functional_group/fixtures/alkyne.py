import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def alkyne(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _alkyne(stk.C(a), stk.C(b), stk.C(c), stk.C(d))


def _alkyne(carbon1, atom1, carbon2, atom2):
    bonders = (carbon1,)
    deleters = (carbon2, atom2)
    return GenericCaseData(
        functional_group=stk.Alkyne(
            carbon1=carbon1,
            atom1=atom1,
            carbon2=carbon2,
            atom2=atom2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon1, atom1, carbon2, atom2),
        bonders=bonders,
        deleters=deleters,
    )
