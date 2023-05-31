import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicSquare(
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
                "FC12C3=C(N=[C+]3)C3=C4C5=C([C+]=N5)[C+]5C6=C7C8=C([C+"
                "]=N8)[C+]3C4(F)C3=C(N=[C+]3)C3=C1C1=C([C+]=N1)[C+]1C("
                "=C(C4=C([C+]=N4)[C+]32)C1(F)C1=C6N=[C+]1)C1=C([C+]=N1"
                ")C75F"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicSquare(
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
                    optimizer=stk.PeriodicCollapser(),
                ),
            ),
            smiles=(
                "FC12C3=C(N=[C+]3)C3=C4C5=C([C+]=N5)[C+]5C6=C7C8=C([C+"
                "]=N8)[C+]3C4(F)C3=C(N=[C+]3)C3=C1C1=C([C+]=N1)[C+]1C("
                "=C(C4=C([C+]=N4)[C+]32)C1(F)C1=C6N=[C+]1)C1=C([C+]=N1"
                ")C75F"
            ),
            name=name,
        ),
    ),
)
def cof_periodic_square(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
