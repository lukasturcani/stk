import os
import stk

test_dir = 'linear_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_construction(amine2, aldehyde2, boronic_acid2, diol2):
    repeat_units = 3
    monomer_joins = 2*repeat_units - 1

    p1 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear(
            repeating_unit='AB',
            orientation=[0, 0],
            n=repeat_units,
            ends='h'
        )
    )

    p2 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear(
            repeating_unit='AB',
            orientation=[1, 1],
            n=repeat_units,
            ends='fg'
        )
    )

    p3 = stk.ConstructedMolecule(
        building_blocks=[boronic_acid2, diol2],
        topology_graph=stk.polymer.Linear(
            repeating_unit='AB',
            orientation=[0, 0],
            n=repeat_units,
            ends='h'
        )
    )

    path = os.path.join(test_dir, 'p1.mol')
    p1.write(path)
    p2.write(path.replace('1', '2'))
    p3.write(path.replace('1', '3'))

    assert p1.bonds_made == monomer_joins
    assert p2.bonds_made == monomer_joins
    assert p3.bonds_made == monomer_joins*2

    monomer_atom_count = len(amine2.atoms) + len(aldehyde2.atoms)

    # 3 atoms are lost at each join in p1 and p2 due to condensation.
    # In p1 an atom is gained due to replacing aldehyde fg with Hs.
    expected_atoms = (
        monomer_atom_count*repeat_units - 3*monomer_joins + 1
    )
    assert len(p1.atoms) == expected_atoms
    assert len(p2.atoms) == expected_atoms - 1

    # p3 loses two atoms because boronic acid OH groups are repalced
    # with H.
    monomer_atom_count = len(boronic_acid2.atoms) + len(diol2.atoms)

    expected_atoms = (
        monomer_atom_count*repeat_units - 6*monomer_joins - 2
    )
    assert len(p3.atoms) == expected_atoms

    assert p1.building_block_counter[amine2] == repeat_units
    assert p1.building_block_counter[aldehyde2] == repeat_units
    assert p2.building_block_counter[amine2] == repeat_units
    assert p2.building_block_counter[aldehyde2] == repeat_units
    assert p3.building_block_counter[boronic_acid2] == repeat_units
    assert p3.building_block_counter[diol2] == repeat_units

    t1 = stk.polymer.Linear('AB', [0, 0], repeat_units, 'h')
    assert repr(p1.topology_graph) == repr(t1)
    t2 = stk.polymer.Linear('AB', [1, 1], repeat_units, 'fg')
    assert repr(p2.topology_graph) == t2
    t3 = stk.polymer.Linear('AB', [0, 0], repeat_units, 'h')
    assert repr(p3.topology_graph) == t3

    expected_bonds = (
        len(amine2.bonds)*repeat_units +
        len(aldehyde2.bonds)*repeat_units -
        monomer_joins*2 + 1
    )
    assert len(p1.bonds) == expected_bonds
    expected_bonds = (
        amine2.mol.GetNumBonds()*repeat_units +
        aldehyde2.mol.GetNumBonds()*repeat_units -
        monomer_joins*2
    )
    assert len(p2.bonds) == expected_bonds

    expected_bonds = (
        len(boronic_acid2.bonds)*repeat_units +
        len(diol2.bonds)*repeat_units -
        monomer_joins*4 - 2
    )
    assert len(p3.bonds) == expected_bonds
