import stk
import pytest


@pytest.fixture(scope='function')
def tmp_rotaxane(tmp_polymer, tmp_macrocycle):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_polymer, tmp_macrocycle],
        topology_graph=stk.rotaxane.NRotaxane('A', 5)
    )


@pytest.fixture(scope='function')
def tmp_rotaxane_alt1(
    tmp_polymer,
    tmp_macrocycle,
    tmp_macrocycle_alt1
):
    return stk.ConstructedMolecule(
        building_blocks=[
            tmp_polymer,
            tmp_macrocycle,
            tmp_macrocycle_alt1
        ],
        topology_graph=stk.rotaxane.NRotaxane('AAB', 5)
    )
