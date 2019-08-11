import pytest
import stk


@pytest.fixture('session')
def six_plus_eight(amine3, aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[amine3, aldehyde4],
        topology_graph=stk.cage.SixPlusEight()
    )


@pytest.fixture
def tmp_six_plus_eight(tmp_amine3, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine3, tmp_aldehyde4],
        topology_graph=stk.cage.SixPlusEight()
    )


@pytest.fixture('session')
def one_plus_one(amine3, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine3, aldehyde3],
        topology_graph=stk.cage.OnePlusOne()
    )


@pytest.fixture
def tmp_one_plus_one(tmp_amine3, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine3, tmp_aldehyde3],
        topology_graph=stk.cage.OnePlusOne()
    )


@pytest.fixture('session')
def two_plus_two(amine3, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine3, aldehyde3],
        topology_graph=stk.cage.TwoPlusTwo()
    )


@pytest.fixture
def tmp_two_plus_two(tmp_amine3, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine3, tmp_aldehyde3],
        topology_graph=stk.cage.TwoPlusTwo()
    )


@pytest.fixture('session')
def four_plus_four(amine3, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine3, aldehyde3],
        topology_graph=stk.cage.FourPlusFour()
    )


@pytest.fixture
def tmp_four_plus_four(tmp_amine3, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine3, tmp_aldehyde3],
        topology_graph=stk.cage.FourPlusFour()
    )


@pytest.fixture('session')
def twelve_plus_thirty(amine2, aldehyde5):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde5],
        topology_graph=stk.cage.TwelvePlusThirty()
    )


@pytest.fixture
def tmp_twelve_plus_thirty(tmp_amine2, tmp_aldehyde5):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde5],
        topology_graph=stk.cage.TwelvePlusThirty()
    )


@pytest.fixture('session')
def two_plus_four(amine2, aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde4],
        topology_graph=stk.cage.TwoPlusFour()
    )


@pytest.fixture
def tmp_two_plus_four(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.TwoPlusFour()
    )


@pytest.fixture('session')
def three_plus_six(amine2, aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde4],
        topology_graph=stk.cage.ThreePlusSix()
    )


@pytest.fixture
def tmp_three_plus_six(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.ThreePlusSix()
    )


@pytest.fixture('session')
def four_plus_eight(amine2, aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde4],
        topology_graph=stk.cage.FourPlusEight()
    )


@pytest.fixture
def tmp_four_plus_eight(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.FourPlusEight()
    )


@pytest.fixture('session')
def five_plus_ten(amine2, aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde4],
        topology_graph=stk.cage.FivePlusTen()
    )


@pytest.fixture
def tmp_five_plus_ten(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.FivePlusTen()
    )


@pytest.fixture('session')
def six_plus_twelve(amine2, aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde4],
        topology_graph=stk.cage.SixPlusTwelve()
    )


@pytest.fixture
def tmp_six_plus_twelve(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.SixPlusTwelve()
    )


@pytest.fixture('session')
def eight_plus_sixteen(amine2, aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde4],
        topology_graph=stk.cage.EightPlusSixteen()
    )


@pytest.fixture
def tmp_eight_plus_sixteen(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.EightPlusSixteen()
    )


@pytest.fixture('session')
def ten_plus_twenty(amine2, aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde4],
        topology_graph=stk.cage.TenPlusTwenty()
    )


@pytest.fixture
def tmp_ten_plus_twenty(tmp_amine2, tmp_aldehyde4):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde4],
        topology_graph=stk.cage.TenPlusTwenty()
    )


@pytest.fixture('session')
def two_plus_three(amine2, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=stk.cage.TwoPlusThree()
    )


@pytest.fixture
def tmp_two_plus_three(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.TwoPlusThree()
    )


@pytest.fixture('session')
def four_plus_six(amine2, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=stk.cage.FourPlusSix()
    )


@pytest.fixture
def tmp_four_plus_six(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.FourPlusSix()
    )


@pytest.fixture('session')
def four_plus_six2(amine2, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=stk.cage.FourPlusSix2()
    )


@pytest.fixture
def tmp_four_plus_six2(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.FourPlusSix2()
    )


@pytest.fixture('session')
def six_plus_nine(amine2, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=stk.cage.SixPlusNine()
    )


@pytest.fixture
def tmp_six_plus_nine(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.SixPlusNine()
    )


@pytest.fixture('session')
def eight_plus_twelve(amine2, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=stk.cage.EightPlusTwelve()
    )


@pytest.fixture
def tmp_eight_plus_twelve(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.EightPlusTwelve()
    )


@pytest.fixture('session')
def twenty_plus_thirty(amine2, aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=stk.cage.TwentyPlusThirty()
    )


@pytest.fixture
def tmp_twenty_plus_thirty(tmp_amine2, tmp_aldehyde3):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde3],
        topology_graph=stk.cage.TwentyPlusThirty()
    )
