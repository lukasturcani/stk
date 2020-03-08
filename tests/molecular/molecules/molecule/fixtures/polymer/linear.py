import stk
import pytest

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrC1=C(Br)[C+]=N1',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='Br[C+]=NC#CBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='AB',
                    num_repeating_units=2
                ),
            ),
            smiles='BrC#CN=[C+]C1=C(C#CN=[C+]C2=C(Br)N=[C+]2)N=[C+]1',
        ),
    ),
)
def polymer_linear(request):
    return request.param
