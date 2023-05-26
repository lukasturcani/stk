import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.Square(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)C(F)(Br)[C+]1Br",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    lattice_size=(2, 2, 1),
                ),
            ),
            smiles=(
                "FC1(Br)C2=C(C3=C([C+]=N3)C3(F)C(=C(C4=C(Br)[C+]=N4)[C"
                "+]3Br)C3=C([C+]=N3)[C+]3C(C4=C(Br)[C+]=N4)=C(C4=C(Br)"
                "[C+]=N4)C3(F)C3=C(N=[C+]3)C3=C(C4=C(Br)[C+]=N4)C(F)(B"
                "r)[C+]3C3=C2N=[C+]3)[C+]1Br"
            ),
            name=name,
        ),
    ),
)
def cof_square(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
