import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def boronic_acid(get_atom_ids):
    a, b, c, d, e, f = get_atom_ids(6)
    return _boronic_acid(
        boron=stk.B(a),
        oxygen1=stk.O(b),
        hydrogen1=stk.H(c),
        oxygen2=stk.O(d),
        hydrogen2=stk.H(e),
        atom=stk.C(f),
    )


def _boronic_acid(boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom):
    bonders = (oxygen1, oxygen2)
    deleters = (hydrogen1, hydrogen2)
    return GenericCaseData(
        functional_group=stk.BoronicAcid(
            boron=boron,
            oxygen1=oxygen1,
            hydrogen1=hydrogen1,
            oxygen2=oxygen2,
            hydrogen2=hydrogen2,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom),
        bonders=bonders,
        deleters=deleters,
    )
