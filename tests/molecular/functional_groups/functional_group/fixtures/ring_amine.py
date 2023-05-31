import pytest
import stk

from ..case_data import CaseData


@pytest.fixture
def ring_amine(get_atom_ids):
    a, b, c, d, e, f, g = get_atom_ids(7)
    return _ring_amine(
        nitrogen=stk.N(a),
        hydrogen1=stk.H(b),
        hydrogen2=stk.H(c),
        carbon1=stk.C(d),
        carbon2=stk.C(e),
        hydrogen3=stk.H(f),
        carbon3=stk.C(g),
    )


def _ring_amine(
    nitrogen,
    hydrogen1,
    hydrogen2,
    carbon1,
    carbon2,
    hydrogen3,
    carbon3,
):
    return CaseData(
        functional_group=stk.RingAmine(
            nitrogen=nitrogen,
            hydrogen1=hydrogen1,
            hydrogen2=hydrogen2,
            carbon1=carbon1,
            carbon2=carbon2,
            hydrogen3=hydrogen3,
            carbon3=carbon3,
        ),
        atoms=(
            nitrogen,
            hydrogen1,
            hydrogen2,
            carbon1,
            carbon2,
            hydrogen3,
            carbon3,
        ),
        placers=(nitrogen,),
        core_atoms=(nitrogen,),
    )
