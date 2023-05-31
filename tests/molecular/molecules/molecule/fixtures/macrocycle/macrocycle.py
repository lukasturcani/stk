import pytest
import stk
from pytest_lazyfixture import lazy_fixture

from ...case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=[C+]1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=4,
                ),
            ),
            smiles=(
                "[C+]1=[C+]C2=C1C1=C(N=[C+]1)C1=C([C+]=[C+]1)C1=C(N=[C"
                "+]1)C1=C([C+]=[C+]1)C1=C(N=[C+]1)C1=C([C+]=[C+]1)C1=C"
                "2N=[C+]1"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
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
            ),
            smiles=(
                "N1=NN2N1N1N=NN1N1N=NN1N1N=NN1N1N=NN1N1N=NN1N1N=NN1N1N" "=NN21"
            ),
            name=name,
        ),
    ),
)
def macrocycle_macrocycle(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )


@pytest.fixture(
    params=(lazy_fixture("macrocycle_macrocycle"),),
)
def macrocycle(request):
    return request.param
