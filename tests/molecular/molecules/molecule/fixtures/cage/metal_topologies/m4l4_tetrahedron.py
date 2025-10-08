import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import get_iron_complex, get_tritopic_linker


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M4L4Tetrahedron(
                    building_blocks={
                        get_iron_complex(): range(4),
                        get_tritopic_linker(): range(4, 8),
                    },
                ),
            ),
            smiles=(
                "[H]C1=C([H])C([H])=[N]2->[Fe+2]3456<-[N]7=C([H])C([H])=C([H]"
                ")C([H])=C7C([H])=[N]->3C3=C([H])C([H])=C(N7C8=C([H])C([H])=C"
                "([N]9=C([H])C%10=C([H])C([H])=C([H])C([H])=[N]%10->[Fe+2]<-9"
                "%10%11%12<-[N]9=C([H])C([H])=C([H])C([H])=C9C([H])=[N]->%10C"
                "9=C([H])C([H])=C(N(C%10=C([H])C([H])=C([N]->4=C([H])C2=C1[H]"
                ")C([H])=C%10[H])C1=C([H])C([H])=C([N]2=C([H])C4=C([H])C([H])"
                "=C([H])C([H])=[N]4->[Fe+2]<-24%10%13<-[N]2=C([H])C([H])=C([H"
                "])C([H])=C2C([H])=[N]->4C2=C([H])C([H])=C(N(C4=C([H])C([H])="
                "C([N]->5=C([H])C5=C([H])C([H])=C([H])C([H])=[N]->65)C([H])=C"
                "4[H])C4=C([H])C([H])=C([N]5=C([H])C6=C([H])C([H])=C([H])C([H"
                "])=[N]6->[Fe+2]<-56%14(<-[N]5=C([H])C([H])=C([H])C([H])=C5C("
                "[H])=[N]->6C5=C([H])C([H])=C7C([H])=C5[H])<-[N]5=C([H])C([H]"
                ")=C([H])C([H])=C5C([H])=[N]->%14C5=C([H])C([H])=C(N(C6=C([H]"
                ")C([H])=C([N]->%11=C([H])C7=C([H])C([H])=C([H])C([H])=[N]->%"
                "127)C([H])=C6[H])C6=C([H])C([H])=C([N]->%10=C([H])C7=C([H])C"
                "([H])=C([H])C([H])=[N]->%137)C([H])=C6[H])C([H])=C5[H])C([H]"
                ")=C4[H])C([H])=C2[H])C([H])=C1[H])C([H])=C9[H])C([H])=C8[H])"
                "C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m4l4_tetrahedron(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
