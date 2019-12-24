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

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('CCCC'),
            factory=stk.AmineFactory(),
            functional_groups=(),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('CNCCNCC'),
            factory=stk.SecondaryAmineFactory(),
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(1), stk.C(0), stk.C(2), stk.H(10)),
                    bonders=(stk.N(1), ),
                    deleters=(stk.H(10), ),
                ),
                stk.Amine(
                    atoms=(stk.N(4), stk.C(3), stk.C(5), stk.H(15)),
                    bonders=(stk.N(4), ),
                    deleters=(stk.H(15), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('O=CCC=O'),
            factory=stk.AldehydeFactory(),
            functional_groups=(
                stk.Aldehyde(
                    atoms=(stk.O(0), stk.C(1), stk.H(5)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.O(0), ),
                ),
                stk.Aldehyde(
                    atoms=(stk.O(4), stk.C(3), stk.H(8)),
                    bonders=(stk.C(3), ),
                    deleters=(stk.O(4), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('O=C(O)CCCC(O)=O'),
            factory=stk.CarboxylicAcidFactory(),
            functional_groups=(
                stk.CarboxylicAcid(
                    atoms=(stk.O(0), stk.C(1), stk.O(2), stk.H(9)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.O(2), stk.H(9)),
                ),
                stk.CarboxylicAcid(
                    atoms=(stk.O(7), stk.O(8), stk.C(6), stk.H(16)),
                    bonders=(stk.C(6), ),
                    deleters=(stk.O(7), stk.H(16)),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('O=C(N)CCC(=O)N'),
            factory=stk.AmideFactory(),
            functional_groups=(
                stk.Amide(
                    atoms=(
                        stk.O(0),
                        stk.C(1),
                        stk.N(2),
                        stk.H(8),
                        stk.H(9),
                    ),
                    bonders=(stk.C(1), ),
                    deleters=(stk.N(2), stk.H(8), stk.H(9)),
                ),
                stk.Amide(
                    atoms=(
                        stk.C(5),
                        stk.O(6),
                        stk.N(7),
                        stk.H(14),
                        stk.H(15),
                    ),
                    bonders=(stk.C(5), ),
                    deleters=(stk.N(7), stk.H(14), stk.H(15)),
                ),
            ),
        ),

    ],
)
def get_functional_groups_test_case(request):
    return request.param
