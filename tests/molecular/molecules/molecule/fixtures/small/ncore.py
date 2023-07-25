import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.small.NCore(
                    core_building_block=(
                        stk.BuildingBlock(
                            smiles="BrC(Br)Br",
                            functional_groups=stk.BromoFactory(),
                        )
                    ),
                    arm_building_blocks=stk.BuildingBlock(
                        smiles="BrC",
                        functional_groups=stk.BromoFactory(),
                    ),
                    repeating_unit="A",
                ),
            ),
            smiles="[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]",
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.small.NCore(
                    core_building_block=(
                        stk.BuildingBlock(
                            smiles="C(Br)1C(Br)C(Br)C(Br)C(Br)C(Br)C1Br",
                            functional_groups=stk.BromoFactory(),
                        )
                    ),
                    arm_building_blocks=stk.BuildingBlock(
                        smiles="BrC",
                        functional_groups=stk.BromoFactory(),
                    ),
                    repeating_unit="A",
                ),
            ),
            smiles=(
                "[H]C([H])([H])C1([H])C([H])(C([H])([H])[H])C([H])(C([H"
                "])([H])[H])C([H])(C([H])([H])[H])C([H])(C([H])([H])[H]"
                ")C([H])(C([H])([H])[H])C1([H])C([H])([H])[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.small.NCore(
                    core_building_block=(
                        stk.BuildingBlock(
                            smiles="C(Br)1C(Br)C(Br)C(Br)C(Br)C(Br)C1Br",
                            functional_groups=stk.BromoFactory(),
                        )
                    ),
                    arm_building_blocks=[
                        stk.BuildingBlock(
                            smiles="BrC",
                            functional_groups=stk.BromoFactory(),
                        ),
                        stk.BuildingBlock(
                            smiles="BrCN",
                            functional_groups=stk.BromoFactory(),
                        ),
                    ],
                    repeating_unit="ABABABA",
                ),
            ),
            smiles=(
                "[H]N([H])C([H])([H])C1([H])C([H])(C([H])([H])[H])C([H]"
                ")(C([H])([H])[H])C([H])(C([H])([H])N([H])[H])C([H])(C("
                "[H])([H])[H])C([H])(C([H])([H])N([H])[H])C1([H])C([H])"
                "([H])[H]"
            ),
            name=name,
        ),
    ),
)
def small_ncore(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
