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
                "[H]C1=C([H])C([H])=N2->[Fe+2]3456<-N7=C([H])C([H])=C"
                "([H])C([H])=C7C([H])=N->3C3=C([H])C([H])=C(C([H])=C3"
                "[H])C3=C([H])C([H])=C(C([H])=C3[H])N3->[Fe+2]789(<-N"
                "%10=C([H])C([H])=C([H])C([H])=C%10C=3[H])<-N3=C([H])"
                "C([H])=C([H])C([H])=C3C([H])=N->7C3=C([H])C([H])=C(C"
                "([H])=C3[H])C3=C([H])C([H])=C(C([H])=C3[H])N3->[Fe+2"
                "]7%10(<-N%11=C([H])C([H])=C([H])C([H])=C%11C([H])=N-"
                ">7C7=C([H])C([H])=C(C([H])=C7[H])C7=C([H])C([H])=C(C"
                "([H])=C7[H])N->4=C([H])C2=C1[H])(<-N1=C([H])C([H])=C"
                "([H])C([H])=C1C=3[H])<-N1=C([H])C([H])=C([H])C([H])="
                "C1C([H])=N->%10C1=C([H])C([H])=C(C([H])=C1[H])C1=C(["
                "H])C([H])=C(C([H])=C1[H])N1->[Fe+2]23(<-N4=C([H])C(["
                "H])=C([H])C([H])=C4C([H])=N->2C2=C([H])C([H])=C(C([H"
                "])=C2[H])C2=C([H])C([H])=C(C([H])=C2[H])N->5=C([H])C"
                "2=C([H])C([H])=C([H])C([H])=N->62)(<-N2=C([H])C([H])"
                "=C([H])C([H])=C2C=1[H])<-N1=C([H])C([H])=C([H])C([H]"
                ")=C1C([H])=N->3C1=C([H])C([H])=C(C([H])=C1[H])C1=C(["
                "H])C([H])=C(C([H])=C1[H])N->8=C([H])C1=C([H])C([H])="
                "C([H])C([H])=N->91"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m4l6_tetrahedron_spacer(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
