def test_is_periodic(case_data):
    _test_is_periodic(case_data.edge, case_data.is_periodic)


def _test_is_periodic(edge, is_periodic):
    assert edge.is_periodic() == is_periodic
