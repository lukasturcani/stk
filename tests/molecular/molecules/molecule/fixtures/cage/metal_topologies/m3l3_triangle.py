import pytest
import stk

from .building_blocks import (
    get_other_linker,
    get_palladium_cispbi_sqpl,
)
from ....case_data import CaseData


@pytest.fixture(
    scope='session',
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L3Triangle(
                    corners=get_palladium_cispbi_sqpl(),
                    linkers=get_other_linker(),
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.GenericFunctionalGroup
                                }): 9,
                            },
                        ),
                    ),
                ),
            ),
            smiles=(
                '[H]C1=C2C([H])=C([H])N(->[Pd+2]3(<-N4=C([H])C([H])=C'
                '(C([H])=C4[H])C4=C([H])C([H])=N(->[Pd+2]5(<-N6=C([H]'
                ')C([H])=C(C([H])=C6[H])C6=C([H])C([H])=N(->[Pd+2]7(<'
                '-N8=C([H])C([H])=C2C([H])=C8[H])<-N([H])([H])C([H])('
                '[H])C([H])([H])N->7([H])[H])C([H])=C6[H])<-N([H])([H'
                '])C([H])([H])C([H])([H])N->5([H])[H])C([H])=C4[H])<-'
                'N([H])([H])C([H])([H])C([H])([H])N->3([H])[H])=C1[H]'
            ),
            name=name,
        ),
    ),
)
def metal_cage_m3l3_triangle(request) -> CaseData:
    return request.param(request.node.originalname)
