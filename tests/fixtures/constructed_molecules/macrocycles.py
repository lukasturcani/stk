import stk
import pytest


@pytest.fixture(scope='function')
def tmp_macrocycle(tmp_bromine2, tmp_bromine2_alt1):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_bromine2, tmp_bromine2_alt1],
        topology_graph=stk.macrocycle.Macrocycle('AB', 3)
    )


@pytest.fixture(scope='function')
def tmp_macrocycle_alt1(tmp_bromine2, tmp_bromine2_alt1):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_bromine2, tmp_bromine2_alt1],
        topology_graph=stk.macrocycle.Macrocycle('AAB', 3)
    )
