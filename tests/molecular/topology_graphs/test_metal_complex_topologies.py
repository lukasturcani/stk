import os
from os.path import join
import numpy as np
import stk
from rdkit.Chem import AllChem as rdkit

test_dir = 'metal_complex_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _build_metal():
    m = rdkit.MolFromSmiles('[Pd+2]')
    m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
    metal = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=None,
    )
    metal_coord_info = {
        0: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        1: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        2: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        3: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
    }
    metal = stk.assign_metal_fgs(
        building_block=metal,
        coordination_info=metal_coord_info
    )
    return metal


def test_metal_definition():
    metal = _build_metal()
    # Test metal position.
    assert np.all(np.isclose(
        metal.get_position_matrix()[0],
        np.array([[0.0, 0.0, 0.0]]),
        rtol=0
    ))
    metal.set_position_matrix(np.array([[1.0, 1.0, 1.0]]))
    assert np.all(np.isclose(
        metal.get_position_matrix()[0],
        np.array([[1.0, 1.0, 1.0]]),
        rtol=0
    ))

    # Test number of FGs and assignment.
    target_fg = stk.FunctionalGroup(
        atoms=tuple([metal.atoms[0]]),
        bonders=tuple([metal.atoms[0]]),
        deleters=(),
        fg_type=stk.FGType(
            name='metal',
            func_group_smarts='',
            bonder_smarts=[],
            deleter_smarts=([])
        ),
    )

    for fg in metal.func_groups:
        assert fg.atoms == target_fg.atoms
        assert fg.bonders == target_fg.bonders
        assert fg.deleters == target_fg.deleters

    # Test metal element.
    for atom in metal.atoms:
        assert atom.atomic_number == 46


def _test_construction(cycle, num_expected_bbs):
    name = cycle.topology_graph.__class__.__name__
    cycle.write(join(test_dir, f'{name}_{len(num_expected_bbs)}.mol'))

    for bb in cycle.get_building_blocks():
        assert cycle.building_block_counter[bb] == num_expected_bbs[bb]

    num_deleters = sum(
        len(fg.deleters)*cycle.building_block_counter[bb]
        for bb in cycle.get_building_blocks() for fg in bb.func_groups
    )
    num_bb_atoms = sum(
        len(bb.atoms)*cycle.building_block_counter[bb]
        for bb in cycle.get_building_blocks()
    )
    num_bb_bonds = sum(
        len(bb.bonds)*cycle.building_block_counter[bb]
        for bb in cycle.get_building_blocks()
    )
    # Check that the correct number of bonds got made.
    assert (
        len(cycle.construction_bonds) ==
        len(cycle.topology_graph.edges)
    )
    # Check correct total number of atoms.
    assert len(cycle.atoms) == num_bb_atoms - num_deleters
    # Check correct total number of bonds.
    assert (
        len(cycle.bonds) ==
        num_bb_bonds + len(cycle.construction_bonds) - num_deleters
    )


def test_sqpl_bidentate_construction():
    metal = _build_metal()
    ligand1 = stk.BuildingBlock(
        'NCCN',
        functional_groups=['amine_metal']
    )

    # Do construction.
    sqpl = stk.metal_complex.SquarePlanarBidentate()
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal, ligand1],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]]),
            ligand1: sqpl.vertices[1:]
        }
    )
    num_expected_bbs = {
        metal: 1,
        ligand1: 2,
    }

    _test_construction(pdl2_sqpl_complex, num_expected_bbs)


def test_sqpl_monodentate_construction():
    metal = _build_metal()
    ligand1 = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    # Handle multiple functional groups.
    ligand1.func_groups = tuple(i for i in [ligand1.func_groups[0]])
    assert len(ligand1.func_groups) == 1

    # Do construction.
    sqpl = stk.metal_complex.SquarePlanarMonodentate()
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal, ligand1],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]]),
            ligand1: sqpl.vertices[1:]
        }
    )
    num_expected_bbs = {
        metal: 1,
        ligand1: 4,
    }

    _test_construction(pdl2_sqpl_complex, num_expected_bbs)


def test_unsat_sqpl_monodentate_construction():
    metal = _build_metal()
    ligand1 = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    # Handle multiple functional groups.
    ligand1.func_groups = tuple(i for i in [ligand1.func_groups[0]])
    assert len(ligand1.func_groups) == 1

    # Do construction.
    sqpl = stk.metal_complex.SquarePlanarMonodentate(
        unsaturated_vertices=[3, 4]
    )
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal, ligand1],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]]),
            ligand1: sqpl.vertices[1:]
        }
    )
    num_expected_bbs = {
        metal: 1,
        ligand1: 2,
    }

    _test_construction(pdl2_sqpl_complex, num_expected_bbs)
