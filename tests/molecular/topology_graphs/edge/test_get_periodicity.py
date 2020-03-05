def test_get_periodicity(case_data):
    _test_get_periodicity(case_data.edge, case_data.periodicity)


def _test_get_periodicity(edge, periodicity):
    assert edge.get_periodicity() == periodicity
