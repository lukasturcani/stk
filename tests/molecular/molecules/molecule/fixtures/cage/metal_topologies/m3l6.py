import numpy as np
import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L6(
                    building_blocks={
                        get_pd_atom(): range(3),
                        get_linker(): range(3, 9),
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9,
                            },
                        ),
                    ),
                ),
            ),
            smiles=(
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H]"
                ")C([H])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->"
                "[Pd+2]9(<-[N]%10=C([H])C([H])=C(C%11=C([H])C(C%12=C([H])"
                "C([H])=[N](->[Pd+2](<-[N]%13=C([H])C([H])=C(C(=C1[H])C=2"
                "[H])C([H])=C%13[H])(<-[N]1=C([H])C([H])=C(C2=C([H])C([H]"
                ")=C([H])C(C%13=C([H])C([H])=[N]->4C([H])=C%13[H])=C2[H])"
                "C([H])=C1[H])<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H]"
                ")C(C4=C([H])C([H])=[N]->9C([H])=C4[H])=C2[H])C([H])=C1[H"
                "])C([H])=C%12[H])=C([H])C([H])=C%11[H])C([H])=C%10[H])<-"
                "[N]1=C([H])C([H])=C(C2=C([H])C(C4=C([H])C([H])=[N]->5C(["
                "H])=C4[H])=C([H])C([H])=C2[H])C([H])=C1[H])C([H])=C8[H])"
                "=C7[H])C([H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L6(
                    building_blocks={
                        get_pd_atom(): range(3),
                        get_linker(): range(3, 9),
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9,
                            },
                        ),
                    ),
                    vertex_positions={0: np.array([2, 0, 2])},
                ),
            ),
            smiles=(
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H])C(["
                "H])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->[Pd+2]9("
                "<-[N]%10=C([H])C([H])=C(C%11=C([H])C(C%12=C([H])C([H])=[N](-"
                ">[Pd+2](<-[N]%13=C([H])C([H])=C(C(=C1[H])C=2[H])C([H])=C%13["
                "H])(<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%13=C([H]"
                ")C([H])=[N]->4C([H])=C%13[H])=C2[H])C([H])=C1[H])<-[N]1=C([H"
                "])C([H])=C(C2=C([H])C([H])=C([H])C(C4=C([H])C([H])=[N]->9C(["
                "H])=C4[H])=C2[H])C([H])=C1[H])C([H])=C%12[H])=C([H])C([H])=C"
                "%11[H])C([H])=C%10[H])<-[N]1=C([H])C([H])=C(C2=C([H])C(C4=C("
                "[H])C([H])=[N]->5C([H])=C4[H])=C([H])C([H])=C2[H])C([H])=C1["
                "H])C([H])=C8[H])=C7[H])C([H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m3l6(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
