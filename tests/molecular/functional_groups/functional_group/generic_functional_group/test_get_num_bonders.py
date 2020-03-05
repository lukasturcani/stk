def test_get_num_bonders(generic_case_data):
    _test_get_num_bonders(
        functional_group=generic_case_data.functional_group,
        num_bonders=len(generic_case_data.bonders),
    )


def _test_get_num_bonders(functional_group, num_bonders):
    assert functional_group.get_num_bonders() == num_bonders
