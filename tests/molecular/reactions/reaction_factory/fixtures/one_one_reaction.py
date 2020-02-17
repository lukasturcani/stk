import pytest
import stk

from ._test_case import _TestCase
from .utilities import MockEdge, get_deleted_atoms


@pytest.fixture(
    params=(

    ),
)
def functional_group(request):
    return request.param


@pytest.fixture
def functional_group1(functional_group):
    return functional_group.clone()


@pytest.fixture
def functional_group2(functional_group):
    return functional_group.clone()


@pytest.fixture(
    params=(
        (0, 0, 0),
        (-1, 0, 1),
        (-1, -1, -2),
        (1, 2, 3),
    ),
)
def periodicity(request):
    return request.param


@pytest.fixture
def one_one_reaction(
    periodicity,
    functional_group1,
    functional_group2,
):
    return _TestCase(
        factory=stk.GenericReactionFactory(),
        construction_state=None,
        edge=MockEdge(periodicity),
        functional_groups=(
            functional_group1,
            functional_group2,
        ),
        reaction_result=stk.ReactionResult(
            new_atoms=(),
            new_bonds=get_new_bonds(
                functional_group1=functional_group1,
                functional_group2=functional_group2,
                periodicity=periodicity,
            ),
            deleted_atoms=get_deleted_atoms(
                functional_group1=functional_group1,
                functional_group2=functional_group2,
            )
        ),
    )


def get_new_bonds(functional_group1, functional_group2, periodicity):
    yield stk.Bond(
        atom1=next(functional_group1.get_bonders()),
        atom2=next(functional_group2.get_bonders()),
        periodicity=periodicity,
    )
