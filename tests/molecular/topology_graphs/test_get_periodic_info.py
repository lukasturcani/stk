import numpy as np
import itertools


def test_get_periodic_info(unscaled_periodic_case):
    """
    Test the collection of the periodic cell from a `.TopologyGraph`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        A test case. Includes the topology graph and expected cell.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_periodic_info(
        topology_graph=unscaled_periodic_case.topology_graph,
        cell=unscaled_periodic_case.cell,
    )


def _test_get_periodic_info(topology_graph, cell):
    """
    Test that the correct cell is extracted.

    Parameters
    ----------
    topology_graph : :class:`.TopologyGraph`
        The topology graph.

    cell : :class:`tuple` of :class:`np.array` or :class:`NoneType`
        Tuple of cell lattice vectors (shape: (3,)) in Angstrom.
        `None` if topology graph is not periodic.

    Returns
    -------
    None : :class:`NoneType`

    """

    actual_cell = (
        topology_graph.get_periodic_info().get_cell_matrix()
    )
    assert np.all(np.array([
        np.allclose(i, j, atol=1e-4)
        for i, j in itertools.zip_longest(actual_cell, cell)
    ]))


def test_get_periodic_info_2(scaled_periodic_case):
    """
    Test getting of :class:`.PeriodicInfo`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The test case. Includes the topology graph and the expected
        cell.

    Returns
    -------
    None : :class:`NoneType`

    """

    construction_result = (
        scaled_periodic_case.topology_graph.construct()
    )
    actual_cell = (
        construction_result.get_periodic_info().get_cell_matrix()
    )
    for i, j in itertools.zip_longest(
        actual_cell,
        scaled_periodic_case.cell,
    ):
        assert np.allclose(i, j, atol=1e-4)
