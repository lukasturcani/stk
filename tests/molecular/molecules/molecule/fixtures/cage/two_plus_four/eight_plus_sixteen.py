import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.EightPlusSixteen(
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
                "[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C(N=[C+]4)[C+]4[C+2][C+"
                "]5C6=C([C+]=N6)[C+]6[C+2][C+]7C8=C([C+]=N8)C38[C+]=N["
                "C+]8[C+]1C1=C([C+]=N1)C13[C+]=N[C+]1[C+]1C8=C([C+]=N8"
                ")C89[C+]=N[C+]8[C+]8C%10=C([C+]=N%10)[C+]7[C+]7N=[C+]"
                "C67C6=C(N=[C+]6)C67[C+]=N[C+]6[C+]6C%10=C([C+]=N%10)["
                "C+]5[C+]5N=[C+]C45C4=C([C+]=N4)[C+]4[C+]2[C+2][C+]2C5"
                "=C([C+]=N5)[C+]3[C+2][C+]1C1=C(N=[C+]1)[C+]1[C+2][C+]"
                "(C3=C(N=[C+]3)C23[C+]=N[C+]43)[C+]2C3=C([C+]=N3)[C+]6"
                "[C+2][C+]7C3=C([C+]=N3)[C+]8[C+2][C+]9C3=C(N=[C+]3)C1"
                "3[C+]=N[C+]23"
            ),
            name=name,
        ),
    ),
)
def cage_eight_plus_sixteen(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
