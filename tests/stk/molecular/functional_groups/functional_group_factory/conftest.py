import pytest
import stk


class GetFunctionalGroupsTestCase:
    def __init__(self, functional_groups, molecule, factory):
        self.functional_groups = functional_groups
        self.molecule = molecule
        self.factory = factory


@pytest.fixture(
    params=[
        GetFunctionalGroupsTestCase(
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(0), stk.H(4), stk.H(5)),
                    bonders=(stk.N(0), ),
                    deleters=(stk.H(4), stk.H(5)),
                ),
                stk.Amine(
                    atoms=(stk.N(3), stk.H(10), stk.H(11)),
                    bonders=(stk.N(3), ),
                    deleters=(stk.H(10), stk.H(11)),
                ),
            ),
            molecule=stk.BuildingBlock('NCCN'),
            factory=stk.AmineFactory(),
        ),
    ],
)
def get_functional_groups_test_case(request):
    return request.param
