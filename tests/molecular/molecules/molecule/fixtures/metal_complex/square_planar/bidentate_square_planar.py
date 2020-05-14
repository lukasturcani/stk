import pytest
import stk

from ....case_data import CaseData


_palladium_atom = stk.BuildingBlock(
    smiles='[Pd+2]',
    functional_groups=(
        stk.SingleAtom(stk.Pd(0, charge=2))
        for i in range(4)
    ),
    position_matrix=([0, 0, 0], ),
)

_bi_1 = stk.BuildingBlock(
    smiles='NCCN',
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#7]~[#6]',
            bonders=(0, ),
            deleters=(),
        ),
    ]
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.BidentateSquarePlanar(
                    metals={_palladium_atom: 0},
                    ligands={_bi_1: (0, 1)},
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom
                                }): 9
                            }
                        )
                    )
                )
            ),
            smiles=(
                '[H]C1([H])C([H])([H])N([H])([H])->[Pd+2]2(<-N1([H])'
                '[H])<-N([H])([H])C([H])([H])C([H])([H])N->2([H])[H]'
            ),
        ),

        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.BidentateSquarePlanar(
                    metals=_palladium_atom,
                    ligands=_bi_1,
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom
                                }): 9
                            }
                        )
                    )
                )
            ),
            smiles=(
                '[H]C1([H])C([H])([H])N([H])([H])->[Pd+2]2(<-N1([H])'
                '[H])<-N([H])([H])C([H])([H])C([H])([H])N->2([H])[H]'
            ),
        ),
    ),
)
def metal_complex_bidentatesquareplanar(request):
    return request.param
