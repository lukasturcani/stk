import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import (
    get_other_linker,
    get_palladium_bi_1,
    get_pd_atom,
)


def _get_palladium_cispbi_sqpl() -> stk.BuildingBlock:
    palladium_cispbi_sqpl = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
            metals=get_pd_atom(),
            ligands=get_palladium_bi_1(),
        ),
    )
    return stk.BuildingBlock.init_from_molecule(
        molecule=palladium_cispbi_sqpl,
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[Pd]~[#7]",
                bonders=(0,),
                deleters=(),
                placers=(0, 1),
            ),
        ],
    )


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.cage.M4L4Square(
                    corners=_get_palladium_cispbi_sqpl(),
                    linkers=get_other_linker(),
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.GenericFunctionalGroup,
                                    }
                                ): 9,
                            }
                        )
                    ),
                )
            ),
            smiles=(
                "[H]C1=C2C([H])=C([H])N(->[Pd+2]3(<-N4=C([H])C([H])=C"
                "(C([H])=C4[H])C4=C([H])C([H])=N(->[Pd+2]5(<-N6=C([H]"
                ")C([H])=C(C([H])=C6[H])C6=C([H])C([H])=N(->[Pd+2]7(<"
                "-N8=C([H])C([H])=C(C([H])=C8[H])C8=C([H])C([H])=N(->"
                "[Pd+2]9(<-N%10=C([H])C([H])=C2C([H])=C%10[H])<-N([H]"
                ")([H])C([H])([H])C([H])([H])N->9([H])[H])C([H])=C8[H"
                "])<-N([H])([H])C([H])([H])C([H])([H])N->7([H])[H])C("
                "[H])=C6[H])<-N([H])([H])C([H])([H])C([H])([H])N->5(["
                "H])[H])C([H])=C4[H])<-N([H])([H])C([H])([H])C([H])(["
                "H])N->3([H])[H])=C1[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M4L4Square(
                    corners=_get_palladium_cispbi_sqpl(),
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
                            }
                        )
                    ),
                    optimizer=stk.MCHammer(
                        num_steps=150,
                        random_seed=1000,
                    ),
                ),
            ),
            smiles=(
                "[H]C1=C2C([H])=C([H])N(->[Pd+2]3(<-N4=C([H])C([H])=C"
                "(C([H])=C4[H])C4=C([H])C([H])=N(->[Pd+2]5(<-N6=C([H]"
                ")C([H])=C(C([H])=C6[H])C6=C([H])C([H])=N(->[Pd+2]7(<"
                "-N8=C([H])C([H])=C(C([H])=C8[H])C8=C([H])C([H])=N(->"
                "[Pd+2]9(<-N%10=C([H])C([H])=C2C([H])=C%10[H])<-N([H]"
                ")([H])C([H])([H])C([H])([H])N->9([H])[H])C([H])=C8[H"
                "])<-N([H])([H])C([H])([H])C([H])([H])N->7([H])[H])C("
                "[H])=C6[H])<-N([H])([H])C([H])([H])C([H])([H])N->5(["
                "H])[H])C([H])=C4[H])<-N([H])([H])C([H])([H])C([H])(["
                "H])N->3([H])[H])=C1[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m4l4_square(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
