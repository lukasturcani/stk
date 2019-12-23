import pytest
import stk


@pytest.fixture(
    params=[

    ],
)
def make_atom_map(request):
    return request.param


class FunctionalGroupAtoms:
    def __init__(self, atoms, bonders, deleters):
        self.atoms = atoms
        self.bonders = bonders
        self.deleters = deleters

    def clone(self):
        return self.__class__(
            atoms=tuple(a.clone() for a in self.atoms),
            bonders=tuple(a.clone() for a in self.bonders),
            deleters=tuple(a.clone() for a in self.deleters),
        )


@pytest.fixture(
    params=[
        FunctionalGroupAtoms(
            atoms=(),
            bonders=(),
            deleters=(),
        ),
        FunctionalGroupAtoms(
            atoms=(),
            bonders=(),
            deleters=(),
        ),
        FunctionalGroupAtoms(
            atoms=(),
            bonders=(),
            deleters=(),
        ),
    ],
)
def functional_group_atoms(request):
    return request.param.clone()


@pytest.fixture(
    params=[

    ],
)
def make_functional_group(request):
    return request.param


@pytest.fixture
def functional_group(make_functional_group):
    return make_functional_group(
        atoms=(stk.C(0), stk.H(2), stk.N(4), stk.P(32)),
        bonders=(stk.C(0), stk.P(32)),
        deleters=(stk.H(2), ),
    )
