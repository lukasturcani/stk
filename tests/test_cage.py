def test_windows():
    ...


def test_window_difference():
    ...


def test_window_variance():
    ...


def test_cavity_size(cc3):
    assert abs(cc3.cavity_size(conformer=1)-5.622072210870494) < 1e-4
