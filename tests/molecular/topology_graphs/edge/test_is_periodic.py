def test_is_periodic(test_case):
    _test_is_periodic(test_case.edge, test_case.is_periodic)


def _test_is_periodic(edge, is_periodic):
    assert edge.is_periodic() == is_periodic
