import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrC1=C(Br)[C+]=[C+]1',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrC1=C(Br)[C+]=N1',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='AB',
                    num_repeating_units=4,
                ),
            ),
            smiles=(
                '[C+]1=[C+]C2=C1C1=C(N=[C+]1)C1=C([C+]=[C+]1)C1=C(N=[C'
                '+]1)C1=C([C+]=[C+]1)C1=C(N=[C+]1)C1=C([C+]=[C+]1)C1=C'
                '2N=[C+]1'
            ),
        ),
    ),
)
def macrocycle(request):
    return request.param
