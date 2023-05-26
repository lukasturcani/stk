import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Honeycomb(
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
                "[H]C12[C+]=NC1(Br)[C+](F)[C+](Br)[C+2][C+]2C1=C([C+]2"
                "[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5=C(N=[C+]5)[C"
                "+]5[C+2][C+](C6=C(Br)[C+]=N6)[C+](F)C6(N=[C+]C56[H])C"
                "5=C([C+]=N5)C56N=[C+]C5([H])[C+](C5=C([C+]7[C+2][C+]("
                "C8=C(Br)[C+]=N8)[C+](F)C8(C9=C(Br)[C+]=N9)N=[C+]C78[H"
                "])N=[C+]5)[C+2][C+](C5=C(N=[C+]5)[C+]5[C+2][C+](C7=C("
                "[C+]=N7)[C+]7[C+2][C+](Br)[C+](F)C8(N=[C+]C78[H])C7=C"
                "(N=[C+]7)C7(N=[C+]C27[H])[C+]3F)C2([H])[C+]=NC2(C2=C("
                "Br)[C+]=N2)[C+]5F)[C+]6F)C2([H])[C+]=NC2(Br)[C+]4F)N="
                "[C+]1"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Honeycomb(
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
                    vertex_alignments={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
                ),
            ),
            smiles=(
                "[H]C12[C+]=NC1(Br)[C+](F)[C+](Br)[C+2][C+]2C1=C([C+]2"
                "[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5=C(N=[C+]5)[C"
                "+]5[C+2][C+](C6=C(Br)[C+]=N6)[C+](F)C6(N=[C+]C56[H])C"
                "5=C([C+]=N5)[C+]5[C+2][C+](C6=C([C+]7[C+2][C+](C8=C(B"
                "r)N=[C+]8)C8([H])[C+]=NC8(C8=C(Br)N=[C+]8)[C+]7F)[C+]"
                "=N6)[C+](F)C6(N=[C+]C56[H])C5=C(N=[C+]5)[C+]5[C+2][C+"
                "](C6=C([C+]=N6)[C+]6[C+2][C+](Br)[C+](F)C7(N=[C+]C67["
                "H])C6=C(N=[C+]6)C6(N=[C+]C26[H])[C+]3F)C2([H])[C+]=NC"
                "2(C2=C(Br)[C+]=N2)[C+]5F)C2([H])[C+]=NC2(Br)[C+]4F)N="
                "[C+]1"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Honeycomb(
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
                    vertex_alignments={0: 2, 1: 2},
                ),
            ),
            smiles=(
                "[H]C12[C+]=NC1(Br)[C+](F)[C+](Br)[C+2][C+]2C1=C([C+]2"
                "[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5=C(N=[C+]5)[C"
                "+]5[C+2][C+](C6=C(Br)[C+]=N6)[C+](F)C6(N=[C+]C56[H])C"
                "5=C([C+]=N5)[C+]5[C+2][C+](C6=C(N=[C+]6)[C+]6[C+2][C+"
                "](C7=C([C+]=N7)[C+]7[C+2][C+](Br)[C+](F)C8(N=[C+]C78["
                "H])C7=C(N=[C+]7)C7(N=[C+]C27[H])[C+]3F)C2([H])[C+]=NC"
                "2(C2=C(Br)[C+]=N2)[C+]6F)C2([H])[C+]=NC2(C2=C(C36N=[C"
                "+]C3([H])[C+](C3=C(Br)[C+]=N3)[C+2][C+](C3=C(Br)[C+]="
                "N3)[C+]6F)N=[C+]2)[C+]5F)C2([H])[C+]=NC2(Br)[C+]4F)N="
                "[C+]1"
            ),
            name=name,
        ),
    ),
)
def cof_honeycomb(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
