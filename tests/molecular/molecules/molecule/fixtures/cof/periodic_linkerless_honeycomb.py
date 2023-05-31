import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicLinkerlessHoneycomb(
                    building_blocks=(
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
                "[H]C12[C+]=NC13[C+](F)[C+]1[C+2][C+]2[C+]2[C+2][C+]4["
                "C+]5[C+2][C+]6[C+]7[C+2][C+]1[C+](F)C1(N=[C+]C71[H])C"
                "17N=[C+]C1([H])[C+]1[C+2][C+]([C+]8[C+2][C+]([C+]9[C+"
                "2][C+]([C+]%10[C+2][C+]1C1([H])[C+]=NC1([C+]%10F)C1(N"
                "=[C+]C61[H])[C+]5F)[C+](F)C1(N=[C+]C91[H])C1(N=[C+]C2"
                "1[H])[C+]4F)C1([H])[C+]=NC13[C+]8F)[C+]7F"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicLinkerlessHoneycomb(
                    building_blocks=(
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
                "[H]C12[C+]=NC13[C+](F)[C+]1[C+2][C+]2[C+]2[C+2][C+]4["
                "C+]5[C+2][C+]6[C+]7[C+2][C+]1[C+](F)C1(N=[C+]C71[H])C"
                "17N=[C+]C1([H])[C+]1[C+2][C+]([C+]8[C+2][C+]([C+]9[C+"
                "2][C+]([C+]%10[C+2][C+]1C1([H])[C+]=NC1([C+]%10F)C1(N"
                "=[C+]C61[H])[C+]5F)[C+](F)C1(N=[C+]C91[H])C1(N=[C+]C2"
                "1[H])[C+]4F)C1([H])[C+]=NC13[C+]8F)[C+]7F"
            ),
            name=name,
        ),
    ),
)
def cof_periodic_linkerless_honeycomb(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
