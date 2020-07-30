import pytest
import numpy as np

from stk.molecular.topology_graphs.cof.cof import NotPeriodicError


def test_get_periodic_cell(periodic_case):
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

    _test_get_periodic_cell(
        topology_graph=periodic_case.topology_graph,
        cell=periodic_case.cell,
    )


def _test_get_periodic_cell(topology_graph, cell):
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

    if not topology_graph._periodic:
        assert cell is None
        with pytest.raises(NotPeriodicError):
            test_cell = topology_graph.get_periodic_cell()

    else:
        test_cell = topology_graph.get_periodic_cell()
        assert np.all(np.array([
            np.allclose(i, j, atol=1e-4)
            for i, j in zip(test_cell, cell)
        ]))
