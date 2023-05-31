import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def difluoro(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _difluoro(stk.F(a), stk.C(b), stk.F(c), stk.C(d))


def _difluoro(fluorine1, atom1, fluorine2, atom2):
    bonders = (atom1, atom2)
    deleters = (fluorine1, fluorine2)
    return GenericCaseData(
        functional_group=stk.Difluoro(
            fluorine1=fluorine1,
            atom1=atom1,
            fluorine2=fluorine2,
            atom2=atom2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(fluorine1, atom1, fluorine2, atom2),
        bonders=bonders,
        deleters=deleters,
    )
