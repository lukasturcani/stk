import pytest
import stk

from ._test_case import _GenericTestCase


@pytest.fixture
def fluoro(get_atom_ids):
    a, b = get_atom_ids(2)
    return _fluoro(stk.F(a), stk.C(b))


def _fluoro(fluorine, atom):
    bonders = (atom, )
    deleters = (fluorine, )
    return _GenericTestCase(
        functional_group=stk.Fluoro(
            fluorine=fluorine,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(fluorine, atom),
        bonders=bonders,
        deleters=deleters,
    )
