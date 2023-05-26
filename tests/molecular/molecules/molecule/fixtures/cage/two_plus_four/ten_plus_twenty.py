import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.TenPlusTwenty(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                "Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+]("
                                "Br)[C+2]1"
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                ),
            ),
            smiles=(
                "[C+]1=NC2=C1C13[C+]=N[C+]1[C+]1C4=C([C+]=N4)[C+]4[C+2"
                "][C+]5C6=C(N=[C+]6)C67[C+]=N[C+]6[C+]6C8=C(N=[C+]8)[C"
                "+]8[C+]9[C+2][C+]%10C%11=C(N=[C+]%11)[C+]%11[C+2][C+]"
                "%12C%13=C(N=[C+]%13)[C+]%13[C+2][C+](C%14=C([C+]=N%14"
                ")C%10%14[C+]=N[C+]8%14)[C+]8C%10=C([C+]=N%10)[C+]6[C+"
                "2][C+]7C6=C(N=[C+]6)[C+]6[C+]([C+2][C+]7C%10=C(N=[C+]"
                "%10)[C+]%12[C+]%10N=[C+]C%11%10C%10=C([C+]=N%10)[C+]%"
                "10[C+2][C+]%11C%12=C([C+]=N%12)[C+]%12[C+2][C+]%14C%1"
                "5=C(N=[C+]%15)C%15%16[C+]=N[C+]%15[C+](C%15=C([C+]=N%"
                "15)C4%15[C+]=N[C+]%15[C+]5C4=C(N=[C+]4)[C+]%14[C+]4N="
                "[C+]C%124C4=C9N=[C+]4)[C+]2[C+2][C+]%16C2=C([C+]=N2)C"
                "%112[C+]=N[C+]2[C+]%10C2=C([C+]=N2)[C+]3[C+2][C+]1C1="
                "C(N=[C+]1)C71[C+]=N[C+]61)C1=C(N=[C+]1)C%131[C+]=N[C+"
                "]81"
            ),
            name=name,
        ),
    ),
)
def cage_ten_plus_twenty(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
