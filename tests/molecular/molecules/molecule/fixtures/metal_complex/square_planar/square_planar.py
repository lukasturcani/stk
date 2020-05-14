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

_mo_1 = stk.BuildingBlock(
    smiles='c1cc2c(cn1)CCCCC2',
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
                stk.metal_complex.SquarePlanar(
                    metals={_palladium_atom: 0},
                    ligands={_mo_1: (0, 1, 2, 3)},
                )
            ),
            smiles=(
                '[H]C1=C([H])N(->[Pd+2](<-N2=C([H])C3=C(C([H])=C2[H])'
                'C([H])([H])C([H])([H])C([H])([H])C([H])([H])C3([H])['
                'H])(<-N2=C([H])C3=C(C([H])=C2[H])C([H])([H])C([H])(['
                'H])C([H])([H])C([H])([H])C3([H])[H])<-N2=C([H])C3=C('
                'C([H])=C2[H])C([H])([H])C([H])([H])C([H])([H])C([H])'
                '([H])C3([H])[H])=C([H])C2=C1C([H])([H])C([H])([H])C(['
                'H])([H])C([H])([H])C2([H])[H]'
            ),
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.SquarePlanar(
                    metals=_palladium_atom,
                    ligands=_mo_1,
                )
            ),
            smiles=(
                '[H]C1=C([H])N(->[Pd+2](<-N2=C([H])C3=C(C([H])=C2[H])'
                'C([H])([H])C([H])([H])C([H])([H])C([H])([H])C3([H])['
                'H])(<-N2=C([H])C3=C(C([H])=C2[H])C([H])([H])C([H])(['
                'H])C([H])([H])C([H])([H])C3([H])[H])<-N2=C([H])C3=C('
                'C([H])=C2[H])C([H])([H])C([H])([H])C([H])([H])C([H])'
                '([H])C3([H])[H])=C([H])C2=C1C([H])([H])C([H])([H])C(['
                'H])([H])C([H])([H])C2([H])[H]'
            ),
        ),
    ),
)
def metal_complex_squareplanar(request):
    return request.param
