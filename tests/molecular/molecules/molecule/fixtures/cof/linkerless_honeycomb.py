import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.LinkerlessHoneycomb(
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
                "[H]C12[C+]=NC1(Br)[C+](F)[C+](Br)[C+2][C+]2[C+]1[C+2]"
                "[C+]2[C+]3[C+2][C+]([C+]4[C+2][C+](Br)[C+](F)C5(N=[C+"
                "]C45[H])C45N=[C+]C4([H])[C+]([C+]4[C+2][C+](Br)[C+](F"
                ")C6(Br)N=[C+]C46[H])[C+2][C+]([C+]4[C+2][C+]([C+]6[C+"
                "2][C+](Br)[C+](F)C7(N=[C+]C67[H])C6(N=[C+]C16[H])[C+]"
                "2F)C1([H])[C+]=NC1(Br)[C+]4F)[C+]5F)C1([H])[C+]=NC1(B"
                "r)[C+]3F"
            ),
            name=name,
        ),
    ),
)
def cof_linkerless_honeycomb(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
