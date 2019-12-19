import pytest
import stk


@pytest.fixture(
    params=[
        stk.ConstructedMolecule(
            building_blocks=[
                stk.BuildingBlock('BrCCBr', ['bromine']),
            ],
            topology_graph=stk.polymer.Linear('A', 3),
        ),
    ],
)
def constructed_molecule(request):
    return request.param.clone()
