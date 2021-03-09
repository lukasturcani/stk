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

palladium_bi_1 = stk.BuildingBlock(
    smiles='[H]N([H])C([H])([H])C([H])([H])N([H])[H]',
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#7]~[#6]',
            bonders=(0, ),
            deleters=(),
        ),
    ]
)
palladium_cispbi_sqpl = stk.ConstructedMolecule(
    stk.metal_complex.CisProtectedSquarePlanar(
        metals=metal_atom,
        ligands=palladium_bi_1,
    )
)
palladium_cispbi_sqpl = stk.BuildingBlock.init_from_molecule(
    molecule=palladium_cispbi_sqpl,
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[Pd]~[#7]',
            bonders=(0, ),
            deleters=(),
            placers=(0, 1),
        ),
    ]
)
linker = stk.BuildingBlock(
    smiles=(
        '[H]C1=NC([H])=C([H])C(C2=C([H])C([H])=NC([H])=C2[H])=C1[H]'
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
                stk.cage.M4L4Square(
                    corners=palladium_cispbi_sqpl,
                    linkers=linker,
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.GenericFunctionalGroup
                                }): 9
                            }
                        )
                    )
                )
            ),
            smiles=(
                '[H]C1=C2C([H])=C([H])N(->[Pd+2]3(<-N4=C([H])C([H])=C'
                '(C([H])=C4[H])C4=C([H])C([H])=N(->[Pd+2]5(<-N6=C([H]'
                ')C([H])=C(C([H])=C6[H])C6=C([H])C([H])=N(->[Pd+2]7(<'
                '-N8=C([H])C([H])=C(C([H])=C8[H])C8=C([H])C([H])=N(->'
                '[Pd+2]9(<-N%10=C([H])C([H])=C2C([H])=C%10[H])<-N([H]'
                ')([H])C([H])([H])C([H])([H])N->9([H])[H])C([H])=C8[H'
                '])<-N([H])([H])C([H])([H])C([H])([H])N->7([H])[H])C('
                '[H])=C6[H])<-N([H])([H])C([H])([H])C([H])([H])N->5(['
                'H])[H])C([H])=C4[H])<-N([H])([H])C([H])([H])C([H])(['
                'H])N->3([H])[H])=C1[H]'
            ),
        ),
    ),
)
def metal_cage_m4l4_square(request):
    return request.param


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.cage.M4L4Square(
                    corners=palladium_cispbi_sqpl,
                    linkers=linker,
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.GenericFunctionalGroup
                                }): 9
                            }
                        )
                    ),
                    optimizer=stk.MCHammer(
                        num_steps=150,
                        random_seed=1000,
                    ),
                )
            ),
            smiles=(
                '[H]C1=C2C([H])=C([H])N(->[Pd+2]3(<-N4=C([H])C([H])=C'
                '(C([H])=C4[H])C4=C([H])C([H])=N(->[Pd+2]5(<-N6=C([H]'
                ')C([H])=C(C([H])=C6[H])C6=C([H])C([H])=N(->[Pd+2]7(<'
                '-N8=C([H])C([H])=C(C([H])=C8[H])C8=C([H])C([H])=N(->'
                '[Pd+2]9(<-N%10=C([H])C([H])=C2C([H])=C%10[H])<-N([H]'
                ')([H])C([H])([H])C([H])([H])N->9([H])[H])C([H])=C8[H'
                '])<-N([H])([H])C([H])([H])C([H])([H])N->7([H])[H])C('
                '[H])=C6[H])<-N([H])([H])C([H])([H])C([H])([H])N->5(['
                'H])[H])C([H])=C4[H])<-N([H])([H])C([H])([H])C([H])(['
                'H])N->3([H])[H])=C1[H]'
            ),
        ),
    ),
)
def metal_cage_opt_m4l4_square(request):
    return request.param
