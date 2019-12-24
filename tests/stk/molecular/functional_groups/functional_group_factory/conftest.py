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

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('O=C(S)CC(S)=O'),
            factory=stk.ThioacidFactory(),
            functional_groups=(
                stk.Thioacid(
                    atoms=(stk.O(0), stk.C(1), stk.S(2), stk.H(7)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.S(2), stk.H(7)),
                ),
                stk.Thioacid(
                    atoms=(stk.C(4), stk.S(5), stk.O(6), stk.H(10)),
                    bonders=(stk.C(4), ),
                    deleters=(stk.S(5), stk.H(10)),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('OCCCO'),
            factory=stk.AlcoholFactory(),
            functional_groups=(
                stk.Alcohol(
                    atoms=(stk.O(0), stk.H(5)),
                    bonders=(stk.O(0), ),
                    deleters=(stk.H(5), ),
                ),
                stk.Alcohol(
                    atoms=(stk.O(4), stk.H(12)),
                    bonders=(stk.O(4), ),
                    deleters=(stk.H(12), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('SCCCS'),
            factory=stk.ThiolFactory(),
            functional_groups=(
                stk.Thiol(
                    atoms=(stk.S(0), stk.H(5)),
                    bonders=(stk.S(0), ),
                    deleters=(stk.H(5), ),
                ),
                stk.Thiol(
                    atoms=(stk.S(4), stk.H(12)),
                    bonders=(stk.S(4), ),
                    deleters=(stk.H(12), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('FCC(F)CCF'),
            factory=stk.FluoroFactory(),
            functional_groups=(
                stk.Fluoro(
                    atoms=(stk.F(0), stk.C(1)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.F(0), ),
                ),
                stk.Fluoro(
                    atoms=(stk.F(3), stk.C(2)),
                    bonders=(stk.C(2), ),
                    deleters=(stk.F(3), ),
                ),
                stk.Fluoro(
                    atoms=(stk.F(6), stk.C(5)),
                    bonders=(stk.C(5), ),
                    deleters=(stk.F(6), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('BrCC(Br)CCBr'),
            factory=stk.BromoFactory(),
            functional_groups=(
                stk.Bromo(
                    atoms=(stk.Br(0), stk.C(1)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.Br(0), ),
                ),
                stk.Bromo(
                    atoms=(stk.Br(3), stk.C(2)),
                    bonders=(stk.C(2), ),
                    deleters=(stk.Br(3), ),
                ),
                stk.Bromo(
                    atoms=(stk.Br(6), stk.C(5)),
                    bonders=(stk.C(5), ),
                    deleters=(stk.Br(6), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('ICC(I)CCI'),
            factory=stk.IodoFactory(),
            functional_groups=(
                stk.Iodo(
                    atoms=(stk.I(0), stk.C(1)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.I(0), ),
                ),
                stk.Iodo(
                    atoms=(stk.I(3), stk.C(2)),
                    bonders=(stk.C(2), ),
                    deleters=(stk.I(3), ),
                ),
                stk.Iodo(
                    atoms=(stk.I(6), stk.C(5)),
                    bonders=(stk.C(5), ),
                    deleters=(stk.I(6), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('C#CC#CC'),
            factory=stk.TerminalAlkyneFactory(),
            functional_groups=(
                stk.TerminalAlkyne(
                    atoms=(stk.C(0), stk.C(1), stk.H(5)),
                    bonders=(stk.C(0), ),
                    deleters=(stk.H(5), ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('C#CC#CC'),
            factory=stk.TerminalAlkyneFactory(delete_carbon=True),
            functional_groups=(
                stk.TerminalAlkyne(
                    atoms=(stk.C(0), stk.C(1), stk.H(5)),
                    bonders=(stk.C(0), ),
                    deleters=(stk.H(5), stk.C(0)),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('C=CC=CC'),
            factory=stk.TerminalAlkeneFactory(),
            functional_groups=(
                stk.TerminalAlkene(
                    atoms=(stk.C(0), stk.C(1), stk.H(5), stk.H(6)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.C(0), stk.H(5), stk.H(6)),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('B(O)(O)CCB(O)O'),
            factory=stk.BoronicAcidFactory(),
            functional_groups=(
                stk.BoronicAcid(
                    atoms=(
                        stk.B(0),
                        stk.O(1),
                        stk.O(2),
                        stk.H(8),
                        stk.H(9),
                    ),
                    bonders=(stk.B(0), ),
                    deleters=(stk.O(1), stk.O(2), stk.H(8), stk.H(9)),
                ),
                stk.BoronicAcid(
                    atoms=(
                        stk.B(5),
                        stk.O(6),
                        stk.O(7),
                        stk.H(14),
                        stk.H(15)
                    ),
                    bonders=(stk.B(5), ),
                    deleters=(
                        stk.O(6),
                        stk.O(7),
                        stk.H(14),
                        stk.H(15),
                    ),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('CC(O)C(O)CC'),
            factory=stk.DiolFactory(),
            functional_groups=(
                stk.Diol(
                    atoms=(
                        stk.C(1),
                        stk.O(2),
                        stk.C(3),
                        stk.O(4),
                        stk.H(11),
                        stk.H(13),
                    ),
                    bonders=(stk.O(2), stk.O(4)),
                    deleters=(stk.H(11), stk.H(13)),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('CC(F)C(F)CC'),
            factory=stk.DifluoroFactory(),
            functional_groups=(
                stk.Difluoro(
                    atoms=(stk.C(1), stk.F(2), stk.C(3), stk.F(4)),
                    bonders=(stk.C(1), stk.C(3)),
                    deleters=(stk.F(2), stk.F(4)),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('CC(Br)C(Br)CC'),
            factory=stk.DibromoFactory(),
            functional_groups=(
                stk.Dibromo(
                    atoms=(stk.C(1), stk.Br(2), stk.C(3), stk.Br(4)),
                    bonders=(stk.C(1), stk.C(3)),
                    deleters=(stk.Br(2), stk.Br(4)),
                ),
            ),
        ),

        GetFunctionalGroupsTestCase(
            molecule=stk.BuildingBlock('NCC(Br)c1c(Br)cccc1'),
            factory=stk.RingAmineFactory(),
            functional_groups=(
                stk.RingAmine(
                    atoms=(
                        stk.N(0),
                        stk.C(1),
                        stk.C(2),
                        stk.C(4),
                        stk.H(11),
                        stk.H(12),
                        stk.H(15),
                    ),
                    bonders=(stk.N(0), stk.C(2)),
                    deleters=(stk.H(11), stk.H(12), stk.H(15)),
                ),
            ),
        ),

    ],
)
def get_functional_groups_test_case(request):
    return request.param
