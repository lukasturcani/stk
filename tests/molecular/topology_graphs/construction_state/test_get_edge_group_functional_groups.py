import itertools as it


def test_get_edge_group_functional_groups(test_case):
    _test_get_edge_group_functional_groups(
        construction_state=test_case.construction_state,
        edge_group=test_case.edge_group,
        functional_groups=test_case.functional_groups,
    )


def _test_get_edge_group_functional_groups(
    construction_state,
    edge_group,
    functional_groups,
):
    functional_groups_ = it.zip_longest(
        construction_state.get_edge_group_functional_groups(
            edge_group=edge_group,
        ),
        functional_groups,
    )
    for fg1, fg2 in functional_groups_:
        fg1 is fg2
