import os
import stk
from collections import Counter

if not os.path.exists('constructed_molecule_tests_output'):
    os.mkdir('constructed_molecule_tests_output')


def test_init(amine2, aldehyde2):
    polymer = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=True
    )
    polymer2 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=True
    )
    assert polymer is polymer2

    polymer3 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=False
    )
    assert polymer is not polymer3

    polymer4 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [1, 0.5], 3),
        use_cache=True
    )
    assert polymer is not polymer4


def test_get_identity_key(
    polymer,
    tmp_polymer,
    amine2,
    aldehyde3,
    amine2_alt1
):
    assert polymer is not tmp_polymer
    assert polymer.get_identity_key() == tmp_polymer.get_identity_key()

    four_plus_six = stk.cage.FourPlusSix()
    cage1 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3, amine2_alt1],
        topology_graph=four_plus_six,
        building_block_vertices={
            aldehyde3: four_plus_six.vertices[:4],
            amine2: four_plus_six.vertices[4:6],
            amine2_alt1: four_plus_six.vertices[6:]
        }
    )
    cage2 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3, amine2_alt1],
        topology_graph=four_plus_six,
        building_block_vertices={
            aldehyde3: four_plus_six.vertices[:4],
            amine2: four_plus_six.vertices[4:5],
            amine2_alt1: four_plus_six.vertices[5:]
        }
    )
    assert cage1.get_identity_key() != cage2.get_identity_key()


def test_get_building_blocks(
    amine2,
    amine2_alt1,
    aldehyde2,
    aldehyde2_alt1,
    aldehyde3,
    polymer,
    polymer_alt1,
    four_plus_six
):
    polymer_bbs = Counter(polymer.get_building_blocks())
    assert len(polymer_bbs) == 2
    assert polymer_bbs[amine2] == 1
    assert polymer_bbs[aldehyde2] == 1

    polymer_alt1_bbs = Counter(polymer_alt1.get_building_blocks())
    assert len(polymer_alt1_bbs) == 2
    assert polymer_alt1_bbs[amine2_alt1] == 1
    assert polymer_alt1_bbs[aldehyde2_alt1] == 1

    four_plus_six_bbs = Counter(four_plus_six.get_building_blocks())
    assert len(four_plus_six_bbs) == 2
    assert four_plus_six_bbs[amine2] == 1
    assert four_plus_six_bbs[aldehyde3] == 1
