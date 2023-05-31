import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def carboxylic_acid(get_atom_ids):
    a, b, c, d, e = get_atom_ids(5)
    return _carboxylic_acid(
        carbon=stk.C(a),
        oxygen1=stk.O(b),
        oxygen2=stk.O(c),
        hydrogen=stk.H(d),
        atom=stk.C(e),
    )


def _carboxylic_acid(carbon, oxygen1, oxygen2, hydrogen, atom):
    bonders = (carbon,)
    deleters = (oxygen2, hydrogen)
    return GenericCaseData(
        functional_group=stk.CarboxylicAcid(
            carbon=carbon,
            oxygen1=oxygen1,
            oxygen2=oxygen2,
            hydrogen=hydrogen,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen1, oxygen2, hydrogen, atom),
        bonders=bonders,
        deleters=deleters,
    )
