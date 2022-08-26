import numpy as np


def test_clone(periodic_case):
    """
    Test :meth:`.PeriodicInfo.clone`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    original = periodic_case.periodic_info
    clone = periodic_case.periodic_info.clone()

    assert np.all(
        np.equal(
            np.array(clone.get_cell_matrix()),
            np.array(original.get_cell_matrix()),
        )
    )
    assert np.all(
        np.equal(
            clone.get_vector_1(),
            original.get_vector_1(),
        )
    )
    assert np.all(
        np.equal(
            clone.get_vector_2(),
            original.get_vector_2(),
        )
    )
    assert np.all(
        np.equal(
            clone.get_vector_3(),
            original.get_vector_3(),
        )
    )
    assert clone.get_a() == original.get_a()
    assert clone.get_b() == original.get_b()
    assert clone.get_c() == original.get_c()
    assert clone.get_alpha() == original.get_alpha()
    assert clone.get_beta() == original.get_beta()
    assert clone.get_gamma() == original.get_gamma()
