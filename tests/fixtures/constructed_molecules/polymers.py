import stk
import pytest


@pytest.fixture(scope='session')
def polymer(amine2, aldehyde2, ab_chain3):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=ab_chain3
    )


@pytest.fixture(scope='function')
def tmp_polymer(tmp_amine2, tmp_aldehyde2, tmp_ab_chain3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde2],
        topology_graph=tmp_ab_chain3
    )


@pytest.fixture('session')
def polymer_alt1(amine2_alt1, aldehyde2_alt1, ab_chain6):
    return stk.ConstructedMolecule(
        building_blocks=[amine2_alt1, aldehyde2_alt1],
        topology_graph=ab_chain6
    )


@pytest.fixture('session')
def polymer_alt2(amine2_alt2, aldehyde2_alt2, ab_chain3):
    return stk.ConstructedMolecule(
        building_blocks=[aldehyde2_alt2, amine2_alt2],
        topology_graph=ab_chain3
    )
