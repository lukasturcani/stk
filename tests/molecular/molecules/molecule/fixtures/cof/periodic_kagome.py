import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicKagome(
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
                    lattice_size=(2, 2, 1),
                ),
            ),
            smiles=(
                "[C+]1=NC2=C1[C+]1[C+]3[C+2][C+]4C5=C(N=[C+]5)C56[C+]="
                "N[C+]5[C+]5C7=C([C+]=N7)[C+]7[C+]8[C+2][C+]9C%10=C(N="
                "[C+]%10)[C+]%10[C+2][C+]%11C%12=C([C+]=N%12)[C+]%12[C"
                "+]%13[C+2][C+]%14C%15=C(N=[C+]%15)C%15%16[C+]=N[C+]%1"
                "5[C+]%15C%17=C([C+]=N%17)[C+]%17[C+]%18[C+2][C+]%19C%"
                "20=C(N=[C+]%20)[C+]%20[C+2][C+]2[C+]2C%21=C([C+]=N%21"
                ")[C+]%21[C+]([C+2][C+](C%22=C(N=[C+]%22)[C+]%16[C+2]["
                "C+]%15C%15=C([C+]=N%15)[C+]%15[C+]([C+2][C+](C%16=C(N"
                "=[C+]%16)C%10%16[C+]=N[C+]%16[C+]%11C%10=C([C+]=N%10)"
                "[C+]%10[C+]([C+2][C+](C%11=C(N=[C+]%11)[C+]6[C+2][C+]"
                "5C5=C([C+]=N5)[C+]5[C+]([C+2][C+](C6=C(N=[C+]6)C%206["
                "C+]=N[C+]26)C2([C+]=N[C+]52)C2=C%18N=[C+]2)C2=C(N=[C+"
                "]2)C92[C+]=N[C+]72)C2([C+]=N[C+]%102)C2=C%13[C+]=N2)C"
                "2=C([C+]=N2)C42[C+]=N[C+]12)C1([C+]=N[C+]%151)C1=C8N="
                "[C+]1)C1=C(N=[C+]1)C%191[C+]=N[C+]%171)C1([C+]=N[C+]%"
                "211)C1=C3[C+]=N1)C1=C([C+]=N1)C%141[C+]=N[C+]%121"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicKagome(
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
                    lattice_size=(2, 2, 1),
                    optimizer=stk.PeriodicCollapser(),
                ),
            ),
            smiles=(
                "[C+]1=NC2=C1[C+]1[C+]3[C+2][C+]4C5=C(N=[C+]5)C56[C+]="
                "N[C+]5[C+]5C7=C([C+]=N7)[C+]7[C+]8[C+2][C+]9C%10=C(N="
                "[C+]%10)[C+]%10[C+2][C+]%11C%12=C([C+]=N%12)[C+]%12[C"
                "+]%13[C+2][C+]%14C%15=C(N=[C+]%15)C%15%16[C+]=N[C+]%1"
                "5[C+]%15C%17=C([C+]=N%17)[C+]%17[C+]%18[C+2][C+]%19C%"
                "20=C(N=[C+]%20)[C+]%20[C+2][C+]2[C+]2C%21=C([C+]=N%21"
                ")[C+]%21[C+]([C+2][C+](C%22=C(N=[C+]%22)[C+]%16[C+2]["
                "C+]%15C%15=C([C+]=N%15)[C+]%15[C+]([C+2][C+](C%16=C(N"
                "=[C+]%16)C%10%16[C+]=N[C+]%16[C+]%11C%10=C([C+]=N%10)"
                "[C+]%10[C+]([C+2][C+](C%11=C(N=[C+]%11)[C+]6[C+2][C+]"
                "5C5=C([C+]=N5)[C+]5[C+]([C+2][C+](C6=C(N=[C+]6)C%206["
                "C+]=N[C+]26)C2([C+]=N[C+]52)C2=C%18N=[C+]2)C2=C(N=[C+"
                "]2)C92[C+]=N[C+]72)C2([C+]=N[C+]%102)C2=C%13[C+]=N2)C"
                "2=C([C+]=N2)C42[C+]=N[C+]12)C1([C+]=N[C+]%151)C1=C8N="
                "[C+]1)C1=C(N=[C+]1)C%191[C+]=N[C+]%171)C1([C+]=N[C+]%"
                "211)C1=C3[C+]=N1)C1=C([C+]=N1)C%141[C+]=N[C+]%121"
            ),
            name=name,
        ),
    ),
)
def cof_periodic_kagome(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
