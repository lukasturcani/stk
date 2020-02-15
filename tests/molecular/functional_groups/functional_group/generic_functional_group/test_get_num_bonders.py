def test_get_num_bonders(generic_test_case):
    _test_get_num_bonders(
        functional_group=generic_test_case.functional_group,
        num_bonders=len(generic_test_case.bonders),
    )


def _test_get_num_bonders(functional_group, num_bonders):
    assert functional_group.get_num_bonders() == num_bonders
