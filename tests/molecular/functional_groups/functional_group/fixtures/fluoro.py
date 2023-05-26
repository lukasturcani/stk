import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def fluoro(get_atom_ids):
    a, b = get_atom_ids(2)
    return _fluoro(stk.F(a), stk.C(b))


def _fluoro(fluorine, atom):
    bonders = (atom,)
    deleters = (fluorine,)
    return GenericCaseData(
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
