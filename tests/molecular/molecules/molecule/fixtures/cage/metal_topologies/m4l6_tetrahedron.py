import pytest
import stk

from ....case_data import CaseData


def _get_fake_complex() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="BrCc1cc(CBr)cc(CBr)c1",
        functional_groups=[stk.BromoFactory()],
    )


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M4L6Tetrahedron(
                    building_blocks={_get_fake_complex(): range(4)},
                ),
            ),
            smiles=(
                "[H]C1=C2C([H])=C3C([H])=C1C([H])([H])C([H])([H])C1=C("
                "[H])C4=C([H])C(=C1[H])C([H])([H])C([H])([H])C1=C([H])"
                "C(=C([H])C(=C1[H])C([H])([H])C3([H])[H])C([H])([H])C("
                "[H])([H])C1=C([H])C(=C([H])C(=C1[H])C([H])([H])C2([H]"
                ")[H])C([H])([H])C4([H])[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m4l6_tetrahedron(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
