import pytest

import stk

from ...case_data import CaseData

# With 1 reaction pair.
cycle1 = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
        building_blocks=(
            stk.BuildingBlock(
                smiles="[Br]CC[Br]",
                functional_groups=[stk.BromoFactory()],
            ),
            stk.BuildingBlock(
                smiles="[Br]C(CCI)C[Br]",
                functional_groups=[stk.BromoFactory()],
            ),
        ),
        repeating_unit="ABAAA",
        num_repeating_units=1,
    ),
)
axle1 = stk.ConstructedMolecule(
    topology_graph=stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            stk.BuildingBlock("BrCNCBr", [stk.BromoFactory()]),
            stk.BuildingBlock("BrCOCI", [stk.BromoFactory()]),
        ),
        repeating_unit="ABABC",
        num_repeating_units=1,
    )
)
rotaxane1 = stk.ConstructedMolecule(
    topology_graph=stk.rotaxane.NRotaxane(
        axle=stk.BuildingBlock.init_from_molecule(axle1),
        cycles=(stk.BuildingBlock.init_from_molecule(cycle1),),
        repeating_unit="A",
        num_repeating_units=1,
    ),
)

# With 2 reaction pairs.
cycle2 = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
        building_blocks=(
            stk.BuildingBlock(
                smiles="[Br]CC[Br]",
                functional_groups=[stk.BromoFactory()],
            ),
            stk.BuildingBlock(
                smiles="[Br]C(CCI)C[Br]",
                functional_groups=[stk.BromoFactory()],
            ),
        ),
        repeating_unit="ABABA",
        num_repeating_units=1,
    ),
)
axle2 = stk.ConstructedMolecule(
    topology_graph=stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            stk.BuildingBlock("BrCNCBr", [stk.BromoFactory()]),
            stk.BuildingBlock("BrCOCI", [stk.BromoFactory()]),
        ),
        repeating_unit="CABABC",
        num_repeating_units=1,
    )
)
rotaxane2 = stk.ConstructedMolecule(
    topology_graph=stk.rotaxane.NRotaxane(
        axle=stk.BuildingBlock.init_from_molecule(axle2),
        cycles=(stk.BuildingBlock.init_from_molecule(cycle2),),
        repeating_unit="A",
        num_repeating_units=1,
    ),
)


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.small.InternalReaction(
                    building_block=stk.BuildingBlock.init_from_molecule(
                        molecule=rotaxane1,
                        functional_groups=(stk.IodoFactory(),),
                    ),
                    num_reactions=1,
                ),
            ),
            smiles=(
                "[H]N(C([H])([H])C([H])([H])OC([H])([H])C([H])([H])C([H])([H])"
                "C1([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H]"
                ")C([H])([H])C([H])([H])C([H])([H])C1([H])[H])C([H])([H])C([H]"
                ")([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H]"
                ")([H])Br"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.small.InternalReaction(
                    building_block=stk.BuildingBlock.init_from_molecule(
                        molecule=rotaxane2,
                        functional_groups=(stk.IodoFactory(),),
                    ),
                    num_reactions=2,
                ),
            ),
            smiles=(
                "[H]N1C([H])([H])C([H])([H])OC([H])([H])C([H])([H])C([H])([H])"
                "C2([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H]"
                ")C([H])(C([H])([H])C([H])([H])C([H])([H])OC([H])([H])C([H])(["
                "H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])(["
                "H])C1([H])[H])C([H])([H])C([H])([H])C2([H])[H]"
            ),
            name=name,
        ),
    ),
)
def small_internal(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
