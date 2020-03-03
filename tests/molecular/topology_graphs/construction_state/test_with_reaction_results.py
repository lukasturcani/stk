import stk
import numpy as np

from .utilities import is_clone


ReactionResult = stk.molecular.reactions.reaction.ReactionResult


def test_with_reaction_results(construction_state):
    clone = construction_state.clone()
    _test_with_reaction_results(construction_state)
    # Test immutability.
    is_clone(construction_state, clone)


def _test_with_reaction_results(old_state):
    results = get_reaction_results(old_state)
    new_state = old_state.with_reaction_results(results)
    check_num_atoms(old_state, new_state, results)
    check_atoms(new_state, results)
    check_position_matrix(new_state, results)
    check_num_bonds(old_state, new_state, results)
    checK_bonds(new_state, results)


def get_reaction_results(construction_state):
    for edge_id in range(0, construction_state.get_num_edges(), 2):
        yield get_reaction_result(construction_state, edge_id)


def get_reaction_result(construction_state, edge_id):
    edge_group = stk.EdgeGroup(
        edges=(construction_state.get_edge(edge_id), ),
    )
    return ReactionResult(
        new_atoms=(
            (stk.C(-1), np.array([0., 0., 0.])),
            (stk.H(-2), np.array([10., 0., 0.])),
        ),
        new_bonds=(
            stk.Bond(stk.C(-1), stk.H(-2), 1),
        ),
        deleted_atoms=tuple(get_atoms(
            construction_state.get_edge_group_functional_groups(
                edge_group=edge_group,
            ),
        )),
    )


def get_atoms(functional_groups):
    for functional_group in functional_groups:
        yield from functional_group.get_atoms()


def check_num_atoms(old_state, new_state, reaction_results):
    num_old_state = sum(1 for _ in old_state.get_atoms())
    num_new_state = sum(1 for _ in new_state.get_atoms())
    num_deleted = sum(
        len(result.deleted_atoms) for result in reaction_results
    )
    num_added = sum(
        len(result.new_atoms) for result in reaction_results
    )
    assert num_new_state == num_old_state + num_added - num_deleted

    num_infos = sum(1 for _ in new_state.get_atom_infos())
    assert num_infos == num_new_state


def check_position_matrix(new_state, reaction_results):
    num_atoms = sum(1 for _ in new_state.get_atoms())
    assert len(new_state.get_position_matrix()) == num_atoms

    num_new_atoms = sum(
        len(result.new_atoms) for result in reaction_results
    )

    position_matrix = new_state.get_position_matrix()
    for atom_id in get_new_atom_ids(new_state):
        assert np.all(np.equal(
            position_matrix[atom_id],

        ))


def get_new_atom_positions(reaction_results):
    for result in reaction_results:
        for atom, position in result.new_atoms:
            yield position
