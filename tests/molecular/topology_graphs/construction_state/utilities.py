import itertools as it

import numpy as np
import stk


def is_clone(construction_state, clone):
    atoms = it.zip_longest(
        construction_state.get_atoms(),
        clone.get_atoms(),
    )
    for atom1, atom2 in atoms:
        assert atom1 is atom2

    bonds = it.zip_longest(
        construction_state.get_bonds(),
        clone.get_bonds(),
    )
    for bond1, bond2 in bonds:
        assert bond1 is bond2

    atom_infos = it.zip_longest(
        construction_state.get_atom_infos(),
        clone.get_atom_infos(),
    )
    for atom_info1, atom_info2 in atom_infos:
        assert atom_info1 is atom_info2

    bond_infos = it.zip_longest(
        construction_state.get_bond_infos(),
        clone.get_bond_infos(),
    )
    for bond_info1, bond_info2 in bond_infos:
        assert bond_info1 is bond_info2

    assert construction_state.get_num_vertices() == clone.get_num_vertices()

    for vertex_id in range(clone.get_num_vertices()):
        assert construction_state.get_vertex(vertex_id) is clone.get_vertex(
            vertex_id
        )

    assert np.all(
        np.equal(
            construction_state.get_position_matrix(),
            clone.get_position_matrix(),
        )
    )

    assert construction_state.get_num_edges() == clone.get_num_edges()
    for edge_id in range(clone.get_num_edges()):
        assert construction_state.get_edge(edge_id) is clone.get_edge(edge_id)

    edge_group = stk.EdgeGroup(
        edges=map(clone.get_edge, range(clone.get_num_edges())),
    )
    functional_groups = it.zip_longest(
        construction_state.get_edge_group_functional_groups(
            edge_group=edge_group,
        ),
        clone.get_edge_group_functional_groups(edge_group),
    )
    for fg1, fg2 in functional_groups:
        assert fg1 is fg2

    counts1 = construction_state.get_building_block_counts()
    counts2 = clone.get_building_block_counts()
    assert (counts1.keys() | counts2.keys()) == counts1.keys()
    for building_block, count1 in counts1.items():
        assert count1 == counts2[building_block]

    assert all(
        np.all(np.equal(actual_constant, expected_constant))
        for actual_constant, expected_constant in it.zip_longest(
            construction_state.get_lattice_constants(),
            clone.get_lattice_constants(),
        )
    )
