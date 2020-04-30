import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def single_atom(get_atom_ids):
    a, = get_atom_ids(1)
    return _single_atom(stk.N(a))


def _single_atom(atom):
    bonders = (atom, )
    deleters = ()
    return GenericCaseData(
        functional_group=stk.SingleAtom(
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(atom, ),
        bonders=bonders,
        deleters=deleters,
    )
