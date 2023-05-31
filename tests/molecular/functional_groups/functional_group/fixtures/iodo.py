import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def iodo(get_atom_ids):
    a, b = get_atom_ids(2)
    return _iodo(stk.I(a), stk.C(b))


def _iodo(iodine, atom):
    bonders = (atom,)
    deleters = (iodine,)
    return GenericCaseData(
        functional_group=stk.Iodo(
            iodine=iodine,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(iodine, atom),
        bonders=bonders,
        deleters=deleters,
    )
