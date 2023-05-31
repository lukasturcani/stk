import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.OnePlusOne(
                    building_blocks=(
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
                "F[C+]1[C+2]C23[C+2][C+]4[C+2][C+]([C+]12)C12[C+2][C+]"
                "(F)[C+]1[C+]3[C+2][C+]4[C+2]2"
            ),
            name=name,
        ),
    ),
)
def cage_one_plus_one(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
