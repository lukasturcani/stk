import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def bromo(get_atom_ids):
    a, b = get_atom_ids(2)
    return _bromo(stk.Br(a), stk.C(b))


def _bromo(bromine, atom):
    bonders = (atom,)
    deleters = (bromine,)
    return GenericCaseData(
        functional_group=stk.Bromo(
            bromine=bromine,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(bromine, atom),
        bonders=bonders,
        deleters=deleters,
    )
