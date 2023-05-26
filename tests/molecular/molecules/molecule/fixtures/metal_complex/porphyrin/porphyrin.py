import pytest
import stk

from ....case_data import CaseData


def _get_zinc_atom() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="[Zn+2]",
        functional_groups=(
            stk.SingleAtom(stk.Zn(0, charge=2)) for i in range(4)
        ),
        position_matrix=([0, 0, 0],),
    )


def _get_quad_1() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=(
            "Brc1ccc(C2=C3C=CC(=C(c4ccc(Br)cc4)C4=NC(=C(c5ccc(Br)cc5)"
            "C5=CC=C([N]5)C(c5ccc(Br)cc5)=C5C=CC2=N5)C=C4)[N]3)cc1"
        ),
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
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
                stk.metal_complex.Porphyrin(
                    metals={_get_zinc_atom(): 0},
                    ligands={_get_quad_1(): 0},
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
                "[H]C1=C([H])C(C2=C3C([H])=C([H])C4=N3->[Zn+2]35<-N"
                "6=C(C([H])=C([H])C6=C(C6=C([H])C([H])=C(Br)C([H])=C6"
                "[H])C6=C([H])C([H])=C2N->36)C(C2=C([H])C([H])=C(Br)C"
                "([H])=C2[H])=C2C([H])=C([H])C(=C4C3=C([H])C([H])=C(Br"
                ")C([H])=C3[H])N->52)=C([H])C([H])=C1Br"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.Porphyrin(
                    metals=_get_zinc_atom(),
                    ligands=_get_quad_1(),
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
                "[H]C1=C([H])C(C2=C3C([H])=C([H])C4=N3->[Zn+2]35<-N"
                "6=C(C([H])=C([H])C6=C(C6=C([H])C([H])=C(Br)C([H])=C6"
                "[H])C6=C([H])C([H])=C2N->36)C(C2=C([H])C([H])=C(Br)C"
                "([H])=C2[H])=C2C([H])=C([H])C(=C4C3=C([H])C([H])=C(Br"
                ")C([H])=C3[H])N->52)=C([H])C([H])=C1Br"
            ),
            name=name,
        ),
    ),
)
def metal_complex_porphyrin(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
