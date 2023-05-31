import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Kagome(
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
                "BrC1=C([C+]2[C+2][C+]3C4=C([C+]=N4)[C+]4[C+2][C+](C5="
                "C([C+]=N5)C56[C+]=N[C+]5[C+](Br)[C+](Br)[C+2][C+]6C5="
                "C(N=[C+]5)C35[C+]=N[C+]5[C+]2C2=C(Br)[C+]=N2)[C+]2C3="
                "C(N=[C+]3)[C+]3[C+]5[C+2][C+]6C7=C([C+]=N7)[C+]7[C+2]"
                "[C+](C8=C([C+]=N8)C89[C+]=N[C+]8[C+](Br)[C+](Br)[C+2]"
                "[C+]9C8=C(N=[C+]8)C68[C+]=N[C+]38)[C+](Br)[C+]3N=[C+]"
                "C73C3=C([C+]=N3)[C+]3[C+2][C+]6C7=C(N=[C+]7)C78[C+]=N"
                "[C+]7[C+](C7=C([C+]=N7)[C+]7[C+]9[C+2][C+](C%10=C(N=["
                "C+]%10)[C+]%10[C+2][C+](C%11=C(Br)[C+]=N%11)[C+](C%11"
                "=C(Br)[C+]=N%11)[C+]%11N=[C+]C%10%11C%10=C([C+]=N%10)"
                "[C+]%10[C+2][C+](C%11=C(N=[C+]%11)C4%11[C+]=N[C+]2%11"
                ")[C+](C2=C5N=[C+]2)[C+]2N=[C+]C%102C2=C9N=[C+]2)C2(C4"
                "=C(Br)[C+]=N4)[C+]=N[C+]72)[C+](C2=C(Br)[C+]=N2)[C+2]"
                "[C+]8C2=C([C+]=N2)[C+]2[C+2][C+](C4=C([C+]=N4)C64[C+]"
                "=N[C+]4[C+]3Br)[C+](Br)[C+]3N=[C+]C23C2=C(Br)[C+]=N2)"
                "N=[C+]1"
            ),
            name=name,
        ),
    ),
)
def cof_kagome(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
