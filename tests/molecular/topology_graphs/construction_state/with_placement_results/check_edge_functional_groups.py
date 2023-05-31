from collections import defaultdict

import stk


def check_edge_functional_groups(
    old_state,
    new_state,
    building_blocks,
    placement_results,
):
    expected = get_expected_edge_functional_groups(
        old_state=old_state,
        building_blocks=building_blocks,
        placement_results=placement_results,
    )

    for edge_id in range(old_state.get_num_edges()):
        edge_group = stk.EdgeGroup((old_state.get_edge(edge_id),))
        expected_fgs = {get_fg_key(fg): fg for fg in expected[edge_id]}
        new_state_fgs = {
            get_fg_key(fg): fg
            for fg in new_state.get_edge_group_functional_groups(edge_group)
        }
        assert expected_fgs.keys() == new_state_fgs.keys()
        for fg_key in expected_fgs:
            is_equivalent_fg(
                fg1=expected_fgs[fg_key],
                fg2=new_state_fgs[fg_key],
            )


def get_expected_edge_functional_groups(
    old_state,
    building_blocks,
    placement_results,
):
    expected = defaultdict(list)
    for edge_id in range(old_state.get_num_edges()):
        edge_group = stk.EdgeGroup((old_state.get_edge(edge_id),))
        expected[edge_id].extend(
            old_state.get_edge_group_functional_groups(edge_group)
        )
    for edge_id, fg in added_functional_groups(
        old_state=old_state,
        building_blocks=building_blocks,
        placement_results=placement_results,
    ):
        expected[edge_id].append(fg)
    return expected


def is_equivalent_fg(fg1, fg2):
    assert fg1.__class__ is fg2.__class__
    assert tuple(fg1.get_atom_ids()) == tuple(fg2.get_atom_ids())
    assert tuple(fg1.get_placer_ids()) == tuple(fg2.get_placer_ids())


def get_fg_key(fg):
    return tuple(fg.get_atom_ids())


def added_functional_groups(
    old_state,
    building_blocks,
    placement_results,
):
    start_atom_id = sum(1 for _ in old_state.get_atoms())
    for building_block, result in zip(
        building_blocks,
        placement_results,
    ):
        atom_map = {
            atom.get_id(): atom.with_id(id_)
            for id_, atom in enumerate(
                building_block.get_atoms(),
                start_atom_id,
            )
        }
        start_atom_id += building_block.get_num_atoms()
        for fg_id, edge_id in result.functional_group_edges.items():
            fg = next(building_block.get_functional_groups(fg_id))
            yield edge_id, fg.with_atoms(atom_map)
