import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def primary_amino(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _primary_amino(stk.N(a), stk.H(b), stk.H(c), stk.C(d))


def _primary_amino(nitrogen, hydrogen1, hydrogen2, atom):
    bonders = (nitrogen,)
    deleters = (hydrogen1, hydrogen2)
    return GenericCaseData(
        functional_group=stk.PrimaryAmino(
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
