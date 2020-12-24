import numpy as np


def test_get_periodic_info(periodic_case):
    """
    Test collection of periodic cell from `.TopologyGraph`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        A test case. Includes the topology graph and expected cell.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_periodic_info(
        topology_graph=periodic_case.topology_graph,
        cell=periodic_case.cell,
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

    try:
        if not topology_graph._periodic:
            assert cell is None
    except AttributeError:
        # Topology must be a PeriodicClass to access periodic info.
        assert hasattr(topology_graph, '_internal')
        test_cell = (
            topology_graph.get_periodic_info().get_cell_matrix()
        )
        assert np.all(np.array([
            np.allclose(i, j, atol=1e-4)
            for i, j in zip(test_cell, cell)
        ]))
