import pytest
import stk

from ._test_case import _GenericTestCase


@pytest.fixture
def bromo(get_atom_ids):
    a, b = get_atom_ids(2)
    return _bromo(stk.Br(a), stk.C(b))


def _bromo(bromine, atom):
    bonders = (atom, )
    deleters = (bromine, )
    return _GenericTestCase(
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
