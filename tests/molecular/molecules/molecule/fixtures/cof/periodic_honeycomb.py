import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicHoneycomb(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                "Br[C+]1C2[C+]=NC2(Br)[C+](F)[C+](Br)" "[C+2]1"
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    lattice_size=(2, 2, 1),
                ),
            ),
            smiles=(
                "[H]C12[C+]=NC13C1=C(N=[C+]1)C14N=[C+]C1([H])[C+]1[C+2"
                "][C+](C5=C([C+]=N5)[C+]5[C+2][C+]6C7=C(N=[C+]7)[C+]7["
                "C+2][C+]8C9=C([C+]=N9)[C+]9[C+2][C+](C%10=C1N=[C+]%10"
                ")C1([H])[C+]=NC1(C1=C(N=[C+]1)C1%10N=[C+]C1([H])[C+]1"
                "[C+2][C+](C%11=C([C+]=N%11)[C+]%11[C+2][C+](C%12=C(N="
                "[C+]%12)[C+]%12[C+2][C+](C%13=C([C+]=N%13)[C+]([C+2]["
                "C+]2C2=C1N=[C+]2)[C+]3F)[C+](F)C1(N=[C+]C%121[H])C1=C"
                "([C+]=N1)C1(N=[C+]C61[H])[C+]5F)C1([H])[C+]=NC1(C1=C("
                "N=[C+]1)C1(N=[C+]C71[H])[C+]8F)[C+]%11F)[C+]%10F)[C+]"
                "9F)[C+]4F"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicHoneycomb(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                "Br[C+]1C2[C+]=NC2(Br)[C+](F)[C+](Br)" "[C+2]1"
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    lattice_size=(2, 2, 1),
                    optimizer=stk.PeriodicCollapser(),
                ),
            ),
            smiles=(
                "[H]C12[C+]=NC13C1=C(N=[C+]1)C14N=[C+]C1([H])[C+]1[C+2"
                "][C+](C5=C([C+]=N5)[C+]5[C+2][C+]6C7=C(N=[C+]7)[C+]7["
                "C+2][C+]8C9=C([C+]=N9)[C+]9[C+2][C+](C%10=C1N=[C+]%10"
                ")C1([H])[C+]=NC1(C1=C(N=[C+]1)C1%10N=[C+]C1([H])[C+]1"
                "[C+2][C+](C%11=C([C+]=N%11)[C+]%11[C+2][C+](C%12=C(N="
                "[C+]%12)[C+]%12[C+2][C+](C%13=C([C+]=N%13)[C+]([C+2]["
                "C+]2C2=C1N=[C+]2)[C+]3F)[C+](F)C1(N=[C+]C%121[H])C1=C"
                "([C+]=N1)C1(N=[C+]C61[H])[C+]5F)C1([H])[C+]=NC1(C1=C("
                "N=[C+]1)C1(N=[C+]C71[H])[C+]8F)[C+]%11F)[C+]%10F)[C+]"
                "9F)[C+]4F"
            ),
            name=name,
        ),
    ),
)
def cof_periodic_honeycomb(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
