import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.TwelvePlusThirty(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                "Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+]("
                                "Br)[C+]1Br"
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                ),
            ),
            smiles=(
                "[C+]1=NC2=C1[C+]1[C+]3C4=C(N=[C+]4)[C+]4[C+]5C6=C([C+"
                "]=N6)[C+]6[C+]2[C+]2C7=C([C+]=N7)[C+]7[C+]8C9=C([C+]="
                "N9)[C+]1[C+]1C9=C(N=[C+]9)[C+]9[C+]%10C%11=C(N=[C+]%1"
                "1)[C+]8[C+]8N=[C+]C8%11C8=C(N=[C+]8)[C+]8[C+]%12C%13="
                "C([C+]=N%13)[C+]%10[C+]%10C%13=C(N=[C+]%13)[C+]%13[C+"
                "]%14C%15=C([C+]=N%15)[C+]%15[C+](C%16=C(N=[C+]%16)C%1"
                "0%16[C+]=N[C+]9%16)[C+](C9=C([C+]=N9)C19[C+]=N[C+]39)"
                "[C+]1C3=C([C+]=N3)[C+]4[C+]3C4=C(N=[C+]4)[C+]4[C+](C9"
                "=C([C+]=N9)C%159[C+]=N[C+]19)[C+]1C9=C([C+]=N9)C9([C+"
                "]=N[C+]%149)[C+]9C%10=C(N=[C+]%10)[C+]%10[C+]%14C%15="
                "C([C+]=N%15)[C+]([C+]8C8=C([C+]=N8)C8%15[C+]=N[C+]8[C"
                "+](C8=C([C+]=N8)[C+]7%11)[C+]7C8=C([C+]=N8)C8([C+]=N["
                "C+]28)[C+]6C2=C(N=[C+]2)[C+]2[C+](C6=C([C+]=N6)C36[C+"
                "]=N[C+]56)[C+]3C5=C([C+]=N5)C5([C+]=N[C+]45)[C+]1C1=C"
                "(N=[C+]1)C1([C+]=N[C+]%101)[C+]1C4=C([C+]=N4)[C+]3C3("
                "[C+]=N[C+]23)C2=C(N=[C+]2)[C+]7[C+]%15C2=C(N=[C+]2)[C"
                "+]%141)C1([C+]=N[C+]%121)C1=C(N=[C+]1)[C+]%139"
            ),
            name=name,
        ),
    ),
)
def cage_twelve_plus_thirty(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
