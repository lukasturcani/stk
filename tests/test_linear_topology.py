import os
import stk
import numpy as np

test_dir = 'linear_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_linear_vertex(tmp_amine2):
    t = stk.polymer.Linear(
        repeating_unit='AB',
        orientation=[0, 0],
        n=3
    )
    v1, v2 = t.vertices[2:4]

    new_coords = v1.place_building_block(tmp_amine2)
    assert np.allclose(
        a=tmp_amine2.get_centroid(tmp_amine2.get_bonder_ids()),
        b=v1.get_position(),
        atol=1e-6
    )
    assert np.allclose(
        a=new_coords,
        b=tmp_amine2.get_position_matrix(),
        atol=1e-6
    )

    new_coords = v2.place_building_block(tmp_amine2)
    assert np.allclose(
        a=tmp_amine2.get_centroid(tmp_amine2.get_bonder_ids()),
        b=v2.get_position(),
        atol=1e-6
    )
    assert np.allclose(
        a=new_coords,
        b=tmp_amine2.get_position_matrix(),
        atol=1e-6
    )


def test_construction(amine2, aldehyde2, boronic_acid2, diol2):
    repeat_units = 3
    monomer_joins = 2*repeat_units - 1

    p1 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear(
            repeating_unit='AB',
            orientation=[1, 1],
            n=repeat_units
        )
    )

    p2 = stk.ConstructedMolecule(
        building_blocks=[boronic_acid2, diol2],
        topology_graph=stk.polymer.Linear(
            repeating_unit='AB',
            orientation=[0, 0],
            n=repeat_units
        )
    )

    path = os.path.join(test_dir, 'p1.mol')
    p1.write(path)
    p2.write(path.replace('1', '2'))

    assert p1.bonds_made == monomer_joins
    assert p2.bonds_made == monomer_joins*2

    num_monomer_atoms = len(amine2.atoms) + len(aldehyde2.atoms)

    # 3 atoms are lost at each join in p1 due to condensation.
    expected_atoms = num_monomer_atoms*repeat_units - 3*monomer_joins
    assert len(p1.atoms) == expected_atoms

    # 6 atoms are lost at each join due to condensation.
    num_monomer_atoms = len(boronic_acid2.atoms) + len(diol2.atoms)
    expected_atoms = num_monomer_atoms*repeat_units - 6*monomer_joins
    assert len(p2.atoms) == expected_atoms

    assert p1.building_block_counter[amine2] == repeat_units
    assert p1.building_block_counter[aldehyde2] == repeat_units
    assert p2.building_block_counter[boronic_acid2] == repeat_units
    assert p2.building_block_counter[diol2] == repeat_units

    t1 = stk.polymer.Linear('AB', [1, 1], repeat_units)
    assert repr(p1.topology_graph) == repr(t1)
    t2 = stk.polymer.Linear('AB', [0, 0], repeat_units)
    assert repr(p2.topology_graph) == repr(t2)

    expected_bonds = (
        len(amine2.bonds)*repeat_units +
        len(aldehyde2.bonds)*repeat_units -
        monomer_joins*2
    )
    assert len(p1.bonds) == expected_bonds
    expected_bonds = (
        len(boronic_acid2.bonds)*repeat_units +
        len(diol2.bonds)*repeat_units -
        monomer_joins*4
    )
    assert len(p2.bonds) == expected_bonds
