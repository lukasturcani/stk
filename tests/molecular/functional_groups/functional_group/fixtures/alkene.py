import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def alkene(get_atom_ids):
    a, b, c, d, e, f = get_atom_ids(6)
    return _alkene(
        carbon1=stk.C(a),
        atom1=stk.C(b),
        atom2=stk.C(c),
        carbon2=stk.C(d),
        atom3=stk.C(e),
        atom4=stk.C(f),
    )


def _alkene(carbon1, atom1, atom2, carbon2, atom3, atom4):
    bonders = (carbon2,)
    deleters = (carbon1, atom1, atom2)
    return GenericCaseData(
        functional_group=stk.Alkene(
            carbon1=carbon1,
            atom1=atom1,
            atom2=atom2,
            carbon2=carbon2,
            atom3=atom3,
            atom4=atom4,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon1, atom1, atom2, carbon2, atom3, atom4),
        bonders=bonders,
        deleters=deleters,
    )
