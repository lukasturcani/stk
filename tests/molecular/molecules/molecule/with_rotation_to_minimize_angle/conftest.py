import pytest
import stk


@pytest.fixture(
    params=(
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('Brc1ccc(Br)cc1Br', [stk.BromoFactory()]),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='BrCNCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=2,
            ),
        ),
    ),
)
def molecule(request):
    """
    A :class:`.Molecule` instance with at least 3 atoms.

    """

    return request.param.clone()
