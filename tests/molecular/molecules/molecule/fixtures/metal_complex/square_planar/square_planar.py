import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_mo_1, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.SquarePlanar(
                    metals={get_pd_atom(): 0},
                    ligands={get_mo_1(): (0, 1, 2, 3)},
                )
            ),
            smiles=(
                "[H]C1=C([H])N(->[Pd+2](<-N2=C([H])C3=C(C([H])=C2[H])"
                "C([H])([H])C([H])([H])C([H])([H])C([H])([H])C3([H])["
                "H])(<-N2=C([H])C3=C(C([H])=C2[H])C([H])([H])C([H])(["
                "H])C([H])([H])C([H])([H])C3([H])[H])<-N2=C([H])C3=C("
                "C([H])=C2[H])C([H])([H])C([H])([H])C([H])([H])C([H])"
                "([H])C3([H])[H])=C([H])C2=C1C([H])([H])C([H])([H])C(["
                "H])([H])C([H])([H])C2([H])[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.SquarePlanar(
                    metals=get_pd_atom(),
                    ligands=get_mo_1(),
                )
            ),
            smiles=(
                "[H]C1=C([H])N(->[Pd+2](<-N2=C([H])C3=C(C([H])=C2[H])"
                "C([H])([H])C([H])([H])C([H])([H])C([H])([H])C3([H])["
                "H])(<-N2=C([H])C3=C(C([H])=C2[H])C([H])([H])C([H])(["
                "H])C([H])([H])C([H])([H])C3([H])[H])<-N2=C([H])C3=C("
                "C([H])=C2[H])C([H])([H])C([H])([H])C([H])([H])C([H])"
                "([H])C3([H])[H])=C([H])C2=C1C([H])([H])C([H])([H])C(["
                "H])([H])C([H])([H])C2([H])[H]"
            ),
            name=name,
        ),
    ),
)
def metal_complex_square_planar(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
