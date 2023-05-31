import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_pd_atom


def _get_bi_1() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="NCCN",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#7]~[#6]",
                bonders=(0,),
                deleters=(),
            ),
        ],
    )


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.BidentateSquarePlanar(
                    metals={get_pd_atom(): 0},
                    ligands={_get_bi_1(): (0, 1)},
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9
                            }
                        )
                    ),
                )
            ),
            smiles=(
                "[H]C1([H])C([H])([H])N([H])([H])->[Pd+2]2(<-N1([H])"
                "[H])<-N([H])([H])C([H])([H])C([H])([H])N->2([H])[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.BidentateSquarePlanar(
                    metals=get_pd_atom(),
                    ligands=_get_bi_1(),
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9
                            }
                        )
                    ),
                )
            ),
            smiles=(
                "[H]C1([H])C([H])([H])N([H])([H])->[Pd+2]2(<-N1([H])"
                "[H])<-N([H])([H])C([H])([H])C([H])([H])N->2([H])[H]"
            ),
            name=name,
        ),
    ),
)
def metal_complex_bidentate_square_planar(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
