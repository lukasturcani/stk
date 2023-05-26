import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_pd_atom


def _get_copper_atom() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="[Cu+2]",
        functional_groups=(
            stk.SingleAtom(stk.Cu(0, charge=2)) for i in range(4)
        ),
        position_matrix=([0, 0, 0],),
    )


def _get_bi_1() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="O=C(O)c1ccc(Br)cc1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#8]~[#1]",
                bonders=(1,),
                deleters=(2,),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#8X1]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )


def _get_bi_2() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="Nc1ccc(C(=O)O)cc1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#8]~[#1]",
                bonders=(1,),
                deleters=(2,),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#8X1]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.Paddlewheel(
                    metals={
                        _get_copper_atom(): (0,),
                        get_pd_atom(): (1,),
                    },
                    ligands={
                        _get_bi_1(): (2, 3),
                        _get_bi_2(): (0, 1),
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
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
                "[H]C1=C([H])C(C2=O->[Pd+2]34<-O=C(C5=C([H])C([H])=C("
                "Br)C([H])=C5[H])O->[Cu+2](<-O2)(<-OC(C2=C([H])C([H])"
                "=C(N([H])[H])C([H])=C2[H])=O->3)<-OC(C2=C([H])C([H])"
                "=C(N([H])[H])C([H])=C2[H])=O->4)=C([H])C([H])=C1Br"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.Paddlewheel(
                    metals=_get_copper_atom(),
                    ligands=_get_bi_1(),
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
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
                "[H]C1=C([H])C(C2=O->[Cu+2]34<-O=C(C5=C([H])C([H])=C("
                "Br)C([H])=C5[H])O->[Cu+2](<-O2)(<-OC(C2=C([H])C([H])"
                "=C(Br)C([H])=C2[H])=O->3)<-OC(C2=C([H])C([H])=C(Br)C"
                "([H])=C2[H])=O->4)=C([H])C([H])=C1Br"
            ),
            name=name,
        ),
    ),
)
def metal_complex_paddlewheel(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
