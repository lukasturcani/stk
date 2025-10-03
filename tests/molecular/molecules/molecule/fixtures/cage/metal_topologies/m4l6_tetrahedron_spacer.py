import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import get_ditopic_linker, get_iron_complex


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M4L6TetrahedronSpacer(
                    building_blocks={
                        get_iron_complex(): range(4),
                        get_ditopic_linker(): range(4, 10),
                    },
                ),
            ),
            smiles=(
                "[H]C1=C([H])C([H])=[N]2->[Fe+2]3456<-[N]7=C([H])C([H])=C([H]"
                ")C([H])=C7C([H])=[N]->3C3=C([H])C([H])=C(C7=C([H])C([H])=C(["
                "N]8=C([H])C9=C([H])C([H])=C([H])C([H])=[N]9->[Fe+2]<-89%10%1"
                "1<-[N]8=C([H])C([H])=C([H])C([H])=C8C([H])=[N]->9C8=C([H])C("
                "[H])=C(C9=C([H])C([H])=C([N]%12=C([H])C%13=C([H])C([H])=C([H"
                "])C([H])=[N]%13->[Fe+2]<-%12%13%14(<-[N]%12=C([H])C([H])=C(["
                "H])C([H])=C%12C([H])=[N]->%13C%12=C([H])C([H])=C(C%13=C([H])"
                "C([H])=C([N]->4=C([H])C2=C1[H])C([H])=C%13[H])C([H])=C%12[H]"
                ")<-[N]1=C([H])C([H])=C([H])C([H])=C1C([H])=[N]->%14C1=C([H])"
                "C([H])=C(C2=C([H])C([H])=C([N]4=C([H])C%12=C([H])C([H])=C([H"
                "])C([H])=[N]%12->[Fe+2]<-4%12%13(<-[N]4=C([H])C([H])=C([H])C"
                "([H])=C4C([H])=[N]->%12C4=C([H])C([H])=C(C%12=C([H])C([H])=C"
                "([N]->5=C([H])C5=C([H])C([H])=C([H])C([H])=[N]->65)C([H])=C%"
                "12[H])C([H])=C4[H])<-[N]4=C([H])C([H])=C([H])C([H])=C4C([H])"
                "=[N]->%13C4=C([H])C([H])=C(C5=C([H])C([H])=C([N]->%10=C([H])"
                "C6=C([H])C([H])=C([H])C([H])=[N]->%116)C([H])=C5[H])C([H])=C"
                "4[H])C([H])=C2[H])C([H])=C1[H])C([H])=C9[H])C([H])=C8[H])C(["
                "H])=C7[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m4l6_tetrahedron_spacer(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
