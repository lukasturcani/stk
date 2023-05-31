import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def amide(get_atom_ids):
    a, b, c, d, e, f = get_atom_ids(6)
    return _amide(
        carbon=stk.C(a),
        oxygen=stk.O(b),
        nitrogen=stk.N(c),
        hydrogen1=stk.H(d),
        hydrogen2=stk.H(e),
        atom=stk.C(f),
    )


def _amide(carbon, oxygen, nitrogen, hydrogen1, hydrogen2, atom):
    bonders = (carbon,)
    deleters = (nitrogen, hydrogen1, hydrogen2)
    return GenericCaseData(
        functional_group=stk.Amide(
            carbon=carbon,
            oxygen=oxygen,
            nitrogen=nitrogen,
            hydrogen1=hydrogen1,
            hydrogen2=hydrogen2,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen, nitrogen, hydrogen1, hydrogen2, atom),
        bonders=bonders,
        deleters=deleters,
    )
