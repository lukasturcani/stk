import stk
import pytest


@pytest.fixture(scope='function')
def tmp_rotaxane(tmp_polymer, tmp_macrocycle):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_polymer, tmp_macrocycle],
        topology_graph=stk.rotaxane.NRotaxane('A', [0], 5)
    )
