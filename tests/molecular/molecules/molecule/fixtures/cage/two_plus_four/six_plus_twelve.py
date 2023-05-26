import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.SixPlusTwelve(
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
                "[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C(N=[C+]4)[C+]4[C+]5[C+"
                "2][C+]6C7=C(N=[C+]7)[C+]7[C+]8[C+2][C+]9C%10=C(N=[C+]"
                "%10)[C+]%10[C+2][C+]%11C%12=C(N=[C+]%12)[C+]%12[C+2]["
                "C+](C%13=C([C+]=N%13)C9%13[C+]=N[C+]7%13)[C+](C7=C([C"
                "+]=N7)C67[C+]=N[C+]47)[C+]4N=[C+]C%124C4=C([C+]=N4)C3"
                "4[C+]=N[C+]4[C+]1C1=C(N=[C+]1)[C+]%11[C+]1N=[C+]C%101"
                "C1=C(N=[C+]1)[C+]1[C+2][C+]2[C+](C2=C5[C+]=N2)[C+]2N="
                "[C+]C12C1=C8[C+]=N1"
            ),
            name=name,
        ),
    ),
)
def cage_six_plus_twelve(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
