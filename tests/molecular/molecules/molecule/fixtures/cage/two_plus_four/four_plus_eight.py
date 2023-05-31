import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusEight(
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
                "[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C(N=[C+]4)C45[C+]=N[C+]"
                "4[C+]4C6=C([C+]=N6)[C+]3[C+]3N=[C+]C13C1=C(N=[C+]1)[C"
                "+]1[C+]2[C+2][C+]2C3=C(N=[C+]3)[C+]3[C+2][C+](C6=C([C"
                "+]=N6)C26[C+]=N[C+]16)[C+]1C2=C([C+]=N2)[C+]4[C+2][C+"
                "]5C2=C(N=[C+]2)C32[C+]=N[C+]12"
            ),
            name=name,
        ),
    ),
)
def cage_four_plus_eight(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
