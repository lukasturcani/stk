import pytest
import stk

from ....case_data import CaseData


metal_atom = stk.BuildingBlock(
    smiles='[Pd+2]',
    functional_groups=(
        stk.SingleAtom(stk.Pd(0, charge=2))
        for i in range(4)
    ),
    position_matrix=([0, 0, 0], ),
)
linker = stk.BuildingBlock(
    smiles=(
        '[H]C1=NC([H])=C([H])C(C2=C([H])C([H])=C([H])C(C3=C([H])C([H]'
        ')=NC([H])=C3[H])=C2[H])=C1[H]'
    ),
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#6]',
            bonders=(1, ),
            deleters=(),
        ),
    ]
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.cage.M3L6(
                    building_blocks={
                        metal_atom: range(3),
                        linker: range(3, 9)
                    },
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
                '[H]C1=C([H])C2=C([H])C(=C1[H])C1=C([H])C([H])=N(->['
                'Pd+2]34<-N5=C([H])C([H])=C(C([H])=C5[H])C5=C([H])C('
                '[H])=C([H])C(=C5[H])C5=C([H])C([H])=N(->[Pd+2]6(<-N'
                '7=C([H])C([H])=C(C([H])=C7[H])C7=C([H])C([H])=C([H]'
                ')C(=C7[H])C7=C([H])C([H])=N(->[Pd+2](<-N8=C([H])C(['
                'H])=C(C([H])=C8[H])C8=C([H])C(=C([H])C([H])=C8[H])C'
                '8=C([H])C([H])=N->3C([H])=C8[H])(<-N3=C([H])C([H])='
                'C(C([H])=C3[H])C3=C([H])C(=C([H])C([H])=C3[H])C3=C('
                '[H])C([H])=N->4C([H])=C3[H])<-N3=C([H])C([H])=C(C(['
                'H])=C3[H])C3=C([H])C(=C([H])C([H])=C3[H])C3=C([H])C'
                '([H])=N->6C([H])=C3[H])C([H])=C7[H])<-N3=C([H])C([H'
                '])=C2C([H])=C3[H])C([H])=C5[H])C([H])=C1[H]'
            ),
        ),
    ),
)
def mcage_m3l6(request):
    return request.param
