import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.rotaxane.NRotaxane(
                    axle=_get_axle(),
                    cycles=(
                        _get_cycle1(),
                        _get_cycle2(),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=3,
                ),
            ),
            smiles=(
                "BrC1=C(C2=C(C3=C(C4=C(C5=C(C6=C(C7=C(C8=C(C9=C(C%10=C"
                "(C%11=C(C%12=C(C%13=C(C%14=C(C%15=C(Br)N=[C+]%15)N=[C"
                "+]%14)N=[C+]%13)N=[C+]%12)N=[C+]%11)N=[C+]%10)N=[C+]9"
                ")N=[C+]8)N=[C+]7)N=[C+]6)N=[C+]5)N=[C+]4)N=[C+]3)N=[C"
                "+]2)N=[C+]1.N1=NN2N1N1N=NN1N1N=NN1N1N=NN1N1N=NN1N1N=N"
                "N1N1N=NN1N1N=NN21.N1=NN2N1N1N=NN1N1N=NN1N1N=NN1N1N=NN"
                "1N1N=NN1N1N=NN1N1N=NN21.N1=NN2N1N1N=NN1N1N=NN1N1N=NN1"
                "N1N=NN1N1N=NN1N1N=NN1N1N=NN21.[C+]1=[C+][C+]2[C+]1[C+"
                "]1[C+]=[C+][C+]1[C+]1[C+]=[C+][C+]1[C+]1[C+]=[C+][C+]"
                "1[C+]1[C+]=[C+][C+]1[C+]1[C+]=[C+][C+]1[C+]1[C+]=[C+]"
                "[C+]1[C+]1[C+]=[C+][C+]21.[C+]1=[C+][C+]2[C+]1[C+]1[C"
                "+]=[C+][C+]1[C+]1[C+]=[C+][C+]1[C+]1[C+]=[C+][C+]1[C+"
                "]1[C+]=[C+][C+]1[C+]1[C+]=[C+][C+]1[C+]1[C+]=[C+][C+]"
                "1[C+]1[C+]=[C+][C+]21.[C+]1=[C+][C+]2[C+]1[C+]1[C+]=["
                "C+][C+]1[C+]1[C+]=[C+][C+]1[C+]1[C+]=[C+][C+]1[C+]1[C"
                "+]=[C+][C+]1[C+]1[C+]=[C+][C+]1[C+]1[C+]=[C+][C+]1[C+"
                "]1[C+]=[C+][C+]21"
            ),
            name=name,
        ),
    ),
)
def rotaxane_nrotaxane(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )


def _get_axle() -> stk.BuildingBlock:
    molecule = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(
                stk.BuildingBlock(
                    smiles="BrC1=C(Br)[C+]=N1",
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles="Br[C+]=NC#CBr",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            repeating_unit="A",
            num_repeating_units=15,
        ),
    )
    return stk.BuildingBlock.init_from_molecule(molecule)


def _get_cycle1() -> stk.BuildingBlock:
    molecule = stk.ConstructedMolecule(
        topology_graph=stk.macrocycle.Macrocycle(
            building_blocks=(
                stk.BuildingBlock(
                    smiles="BrN1N(Br)N=N1",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            repeating_unit="A",
            num_repeating_units=8,
        ),
    )

    return stk.BuildingBlock.init_from_molecule(molecule)


def _get_cycle2() -> stk.BuildingBlock:
    molecule = stk.ConstructedMolecule(
        topology_graph=stk.macrocycle.Macrocycle(
            building_blocks=(
                stk.BuildingBlock(
                    smiles="Br[C+]1[C+](Br)[C+]=[C+]1",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            repeating_unit="A",
            num_repeating_units=8,
        ),
    )
    return stk.BuildingBlock.init_from_molecule(molecule)
