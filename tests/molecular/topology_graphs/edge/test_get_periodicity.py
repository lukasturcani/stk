def test_get_periodicity(test_case):
    _test_get_periodicity(test_case.edge, test_case.periodicity)


def _test_get_periodicity(edge, periodicity):
    assert edge.get_periodicity() == periodicity
