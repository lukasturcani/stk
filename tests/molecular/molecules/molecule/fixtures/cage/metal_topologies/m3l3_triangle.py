import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import (
    get_other_linker,
    get_palladium_cispbi_sqpl,
)


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L3Triangle(
                    corners=get_palladium_cispbi_sqpl(),
                    linkers=get_other_linker(),
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.GenericFunctionalGroup,
                                    }
                                ): 9,
                            },
                        ),
                    ),
                ),
            ),
            smiles=(
                "[H]C1=C2C([H])=C([H])N(->[Pd+2]3(<-N4=C([H])C([H])=C"
                "(C([H])=C4[H])C4=C([H])C([H])=N(->[Pd+2]5(<-N6=C([H]"
                ")C([H])=C(C([H])=C6[H])C6=C([H])C([H])=N(->[Pd+2]7(<"
                "-N8=C([H])C([H])=C2C([H])=C8[H])<-N([H])([H])C([H])("
                "[H])C([H])([H])N->7([H])[H])C([H])=C6[H])<-N([H])([H"
                "])C([H])([H])C([H])([H])N->5([H])[H])C([H])=C4[H])<-"
                "N([H])([H])C([H])([H])C([H])([H])N->3([H])[H])=C1[H]"
            ),
            name=name,
        ),
        # Non-metal-containing tests with increasing angles between
        # functional groups.
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L3Triangle(
                    corners=stk.BuildingBlock(
                        smiles="C1CC(C1Br)Br",
                        functional_groups=[stk.BromoFactory()],
                    ),
                    linkers=stk.BuildingBlock(
                        smiles="C(#CBr)Br",
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
            ),
            smiles=(
                "[H]C1([H])C([H])([H])C2([H])C#CC3([H])C([H])([H])C([H"
                "])([H])C3([H])C#CC3([H])C([H])([H])C([H])([H])C3([H])"
                "C#CC12[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L3Triangle(
                    corners=stk.BuildingBlock(
                        smiles="C1CC(CC(C1)Br)Br",
                        functional_groups=[stk.BromoFactory()],
                    ),
                    linkers=stk.BuildingBlock(
                        smiles="C(#CBr)Br",
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
            ),
            smiles=(
                "[H]C1([H])C([H])([H])C2([H])C#CC3([H])C([H])([H])C([H"
                "])([H])C([H])([H])C([H])(C#CC4([H])C([H])([H])C([H])("
                "[H])C([H])([H])C([H])(C#CC([H])(C1([H])[H])C2([H])[H]"
                ")C4([H])[H])C3([H])[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L3Triangle(
                    corners=stk.BuildingBlock(
                        smiles="C(CBr)Br",
                        functional_groups=[stk.BromoFactory()],
                    ),
                    linkers=stk.BuildingBlock(
                        smiles="C(#CBr)Br",
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
            ),
            smiles=(
                "[H]C1([H])C#CC([H])([H])C([H])([H])C#CC([H])([H])C([H"
                "])([H])C#CC1([H])[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m3l3_triangle(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
