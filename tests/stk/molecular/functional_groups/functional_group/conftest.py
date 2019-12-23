import pytest
import stk


def make_atom_map_0(functional_group):
    atoms = (stk.Li(200), )
    return dict(zip(functional_group.get_atom_ids(), atoms))


@pytest.fixture(
    params=[
        lambda functional_group: None,
        lambda functional_group: {},
        make_atom_map_0,
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
            atoms=(stk.N(0), stk.H(21), stk.K(3)),
            bonders=(stk.N(0), ),
            deleters=(stk.K(3), ),
        ),
        FunctionalGroupAtoms(
            atoms=(stk.H(120), stk.K(32), stk.H(1)),
            bonders=(stk.K(32), ),
            deleters=(),
        ),
        FunctionalGroupAtoms(
            atoms=(stk.C(7), stk.N(12), stk.Li(2)),
            bonders=(),
            deleters=(stk.N(12), ),
        ),
    ],
)
def functional_group_atoms(request):
    return request.param.clone()


class FunctionalGroupData:
    def __init__(self, functional_group, atoms, bonders, deleters):
        self.functional_group = functional_group
        self.atoms = atoms
        self.bonders = bonders
        self.deleters = deleters


@pytest.fixture(
    params=[
        stk.Amine,
        stk.Aldehyde,
        stk.CarboxylicAcid,
        stk.Amide,
        stk.Thioacid,
        stk.Alcohol,
        stk.Thiol,
        stk.Fluoro,
        stk.Bromo,
        stk.Iodo,
        stk.TerminalAlkyne,
        stk.TerminalAlkene,
        stk.BoronicAcid,
        stk.Diol,
        stk.Difluoro,
        stk.Dibromo,
        stk.RingAmine,
    ],
)
def make_functional_group_0(request, functional_group_atoms):

    def inner():
        return FunctionalGroupData(
            functional_group=request.param(
                atoms=functional_group_atoms.atoms,
                bonders=functional_group_atoms.bonders,
                deleters=functional_group_atoms.deleters,
            ),
            atoms=functional_group_atoms.atoms,
            bonders=functional_group_atoms.bonders,
            deleters=functional_group_atoms.deleters,
        )

    return inner


@pytest.fixture(
    params=[
        pytest.lazy_fixture('make_functional_group_0'),
    ],
)
def make_functional_group(request):
    return request.param


@pytest.fixture
def functional_group(make_functional_group):
    return make_functional_group().functional_group
