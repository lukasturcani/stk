import stk


def get_metal_atom() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles='[Pd+2]',
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2))
            for i in range(4)
        ),
        position_matrix=([0, 0, 0], ),
    )


def get_linker() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=(
            '[H]C1=NC([H])=C([H])C(C2=C([H])C([H])=C([H])C(C3=C([H])'
            'C([H])=NC([H])=C3[H])=C2[H])=C1[H]'
        ),
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts='[#6]~[#7X2]~[#6]',
                bonders=(1, ),
                deleters=(),
            ),
        ]
    )
