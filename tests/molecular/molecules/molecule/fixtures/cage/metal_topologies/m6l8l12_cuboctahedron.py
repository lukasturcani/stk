import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M6L8L12Cube(
                    building_blocks={
                        get_pd_atom(): range(6),
                        get_linker(): range(6, 18),
                        get_linker(): range(18, 26),
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
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
                "[H]C1=[O+][Zr]234567OC([H])=[O+][Zr]89%10%11%12%13O"
                "C([H])=[O+][Zr]%14%15%16%17%18%19OC([H])=[O+][Hf]%2"
                "0%21([O+]=C([H])O[Zr]%22%23(O1)(OC([H])=[O+]%14)(OC"
                "([H])=[O+][Hf]([O+]=C([H])O2)([O+]=C([H])O8)([O+]=C"
                "([H])O%15)([O+]%223)([O+]49)([O+]%23%16)[O+]%10%17)"
                "([O+]5%20)[O+]%18%21)([O+]=C([H])O6)([O+]=C([H])O%1"
                "1)([O+]7%12)[O+]%13%19"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M6L12Cube(
                    building_blocks={
                        get_pd_atom(): range(6),
                        get_linker(): range(6, 18),
                        get_linker(): range(18, 26),
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
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
                    scale_multiplier=1.2,
                ),
            ),
            smiles=(
                "[H]C1=[O+][Zr]234567OC([H])=[O+][Zr]89%10%11%12%13O"
                "C([H])=[O+][Zr]%14%15%16%17%18%19OC([H])=[O+][Hf]%2"
                "0%21([O+]=C([H])O[Zr]%22%23(O1)(OC([H])=[O+]%14)(OC"
                "([H])=[O+][Hf]([O+]=C([H])O2)([O+]=C([H])O8)([O+]=C"
                "([H])O%15)([O+]%223)([O+]49)([O+]%23%16)[O+]%10%17)"
                "([O+]5%20)[O+]%18%21)([O+]=C([H])O6)([O+]=C([H])O%1"
                "1)([O+]7%12)[O+]%13%19"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m6l8l12_cuboctahedron(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
