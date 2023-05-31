import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def secondary_amino(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _secondary_amino(stk.N(a), stk.H(b), stk.C(c), stk.C(d))


def _secondary_amino(nitrogen, hydrogen, atom1, atom2):
    bonders = (nitrogen,)
    deleters = (hydrogen,)
    return GenericCaseData(
        functional_group=stk.SecondaryAmino(
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
