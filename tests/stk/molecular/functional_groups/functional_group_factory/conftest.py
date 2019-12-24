import pytest
import stk


class GetFunctionalGroupsTestCase:
    def __init__(self, molecule, factory, functional_groups):
        self.molecule = molecule
        self.factory = factory
        self.functional_groups = functional_groups


@pytest.fixture(
    params=[
        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('NCCN'),
            factory=stk.AmineFactory(),
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
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('NCCN'),
            factory=stk.AmineFactory(num_deleters=1),
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(0), stk.H(4), stk.H(5)),
                    bonders=(stk.N(0), ),
                    deleters=(stk.H(4), ),
                ),
                stk.Amine(
                    atoms=(stk.N(3), stk.H(10), stk.H(11)),
                    bonders=(stk.N(3), ),
                    deleters=(stk.H(10), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('NCCN'),
            factory=stk.AmineFactory(num_deleters=0),
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(0), stk.H(4), stk.H(5)),
                    bonders=(stk.N(0), ),
                    deleters=(),
                ),
                stk.Amine(
                    atoms=(stk.N(3), stk.H(10), stk.H(11)),
                    bonders=(stk.N(3), ),
                    deleters=(),
                ),
            ),
        ),

    ],
)
def get_functional_groups_test_case(request):
    return request.param
