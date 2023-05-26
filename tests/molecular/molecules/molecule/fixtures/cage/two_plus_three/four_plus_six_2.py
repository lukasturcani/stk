import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix2(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                "Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]"
                                "C2(Br)[C+2]1"
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                ),
            ),
            smiles=(
                "F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)C56[C+"
                "2][C+]7[C+2][C+](C8=C(N=[C+]8)C89[C+2][C+]([C+2][C+]("
                "C%10=C([C+]=N%10)C%10%11[C+2][C+]([C+2][C+](C%12=C2[C"
                "+]=N%12)[C+]%10[C+](F)[C+2]%11)C2=C4[C+]=N2)[C+]8[C+]"
                "(F)[C+2]9)C2=C7[C+]=N2)[C+]5[C+](F)[C+2]6)[C+]13"
            ),
            name=name,
        ),
    ),
)
def cage_four_plus_six_2(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
