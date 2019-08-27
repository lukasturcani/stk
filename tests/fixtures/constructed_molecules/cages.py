import pytest
import stk
from os.path import join

data_dir = join('..', 'data')


@pytest.fixture(scope='function')
def tmp_cc3():
    amine = stk.BuildingBlock.init_from_file(
        path=join(data_dir, 'cc3_amine.mol'),
        functional_groups=['amine']
    )
    aldehyde = stk.BuildingBlock.init_from_file(
        path=join(data_dir, 'cc3_aldehyde.mol'),
        functional_groups=['aldehyde']
    )
    return stk.ConstructedMolecule(
        building_blocks=[amine, aldehyde],
        topology_graph=stk.cage.FourPlusSix()
    )


@pytest.fixture(scope='function')
def tmp_opt_cc3():
    amine = stk.BuildingBlock.init_from_file(
        path=join(data_dir, 'cc3_amine.mol'),
        functional_groups=['amine']
    )
    aldehyde = stk.BuildingBlock.init_from_file(
        path=join(data_dir, 'cc3_aldehyde.mol'),
        functional_groups=['aldehyde']
    )
    cc3 = stk.ConstructedMolecule(
        building_blocks=[amine, aldehyde],
        topology_graph=stk.cage.FourPlusSix()
    )
    return cc3.update_from_file(join(data_dir, 'cc3.mol'))


@pytest.fixture(scope='function')
def tmp_six_plus_eight(tmp_amine3, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine3, tmp_aldehyde4],
        topology_graph=stk.cage.SixPlusEight()
    )


@pytest.fixture(scope='function')
def tmp_one_plus_one(tmp_amine3, tmp_aldehyde3):
    topology_graph = stk.cage.OnePlusOne()
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine3, tmp_aldehyde3],
        topology_graph=topology_graph,
        building_block_vertices={
            tmp_amine3: topology_graph.vertices[:1],
            tmp_aldehyde3: topology_graph.vertices[1:]
        }
    )


@pytest.fixture(scope='function')
def tmp_two_plus_two(tmp_amine3, tmp_aldehyde3):
    topology_graph = stk.cage.TwoPlusTwo()
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine3, tmp_aldehyde3],
        topology_graph=topology_graph,
        building_block_vertices={
            tmp_amine3: topology_graph.vertices[:2],
            tmp_aldehyde3: topology_graph.vertices[2:]
        }
    )


@pytest.fixture(scope='function')
def tmp_four_plus_four(tmp_amine3, tmp_aldehyde3):
    topology_graph = stk.cage.FourPlusFour()
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine3, tmp_aldehyde3],
        topology_graph=topology_graph,
        building_block_vertices={
            tmp_amine3: topology_graph.vertices[:4],
            tmp_aldehyde3: topology_graph.vertices[4:]
        }
    )


@pytest.fixture(scope='function')
def tmp_twelve_plus_thirty(tmp_amine2, tmp_aldehyde5):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde5],
        topology_graph=stk.cage.TwelvePlusThirty()
    )


@pytest.fixture(scope='function')
def tmp_two_plus_four(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.TwoPlusFour()
    )


@pytest.fixture(scope='function')
def tmp_three_plus_six(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.ThreePlusSix()
    )


@pytest.fixture(scope='function')
def tmp_four_plus_eight(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.FourPlusEight()
    )


@pytest.fixture(scope='function')
def tmp_five_plus_ten(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.FivePlusTen()
    )


@pytest.fixture(scope='function')
def tmp_six_plus_twelve(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.SixPlusTwelve()
    )


@pytest.fixture(scope='function')
def tmp_eight_plus_sixteen(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.EightPlusSixteen()
    )


@pytest.fixture(scope='function')
def tmp_ten_plus_twenty(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.TenPlusTwenty()
    )


@pytest.fixture(scope='function')
def tmp_two_plus_three(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.TwoPlusThree()
    )


@pytest.fixture(scope='function')
def four_plus_six(amine2, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=stk.cage.FourPlusSix()
    )


@pytest.fixture(scope='function')
def tmp_four_plus_six(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.FourPlusSix()
    )


@pytest.fixture(scope='function')
def tmp_four_plus_six2(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.FourPlusSix2()
    )


@pytest.fixture(scope='function')
def tmp_six_plus_nine(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.SixPlusNine()
    )


@pytest.fixture(scope='function')
def tmp_eight_plus_twelve(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.EightPlusTwelve()
    )


@pytest.fixture(scope='function')
def tmp_twenty_plus_thirty(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.TwentyPlusThirty()
    )
