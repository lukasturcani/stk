import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def metal_bound_atom(get_atom_ids):
    a, b = get_atom_ids(2)
    return _metal_bound_atom(stk.N(a), stk.Fe(b))


def _metal_bound_atom(atom, metal):
    atoms = (atom, metal)
    binders = (atom, )
    deleters = ()
    return GenericCaseData(
        functional_group=stk.MetalBoundAtom(atom=atom, metal=metal),
        atoms=atoms,
        bonders=binders,
        deleters=deleters,
    )
