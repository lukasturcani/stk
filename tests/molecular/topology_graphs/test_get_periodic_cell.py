import pytest
import numpy as np

from stk.molecular.topology_graphs.cof.cof import NotPeriodicError


def test_get_periodic_cell(periodic_case):
    """
    Test correction deletion of atoms by a :class:`.Reaction`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Includes the reaction and atoms which should have
        been deleted.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_periodic_cell(
        topology_graph=periodic_case.topology_graph,
        periodic_cell=periodic_case.cell,
    )


def _test_get_periodic_cell(topology_graph, periodic_cell):
    """
    Test that the correct atoms are deleted by `reaction_result`.

    Parameters
    ----------
    reaction_result : :class:`.ReactionResult`
        The result of a reaction.

    deleted_atoms : :class:`tuple` of :class:`.Atom`
        The atoms which should be deleted by the reaction.

    Returns
    -------
    None : :class:`NoneType`

    """

    if not topology_graph._periodic:
        assert periodic_cell is None
        with pytest.raises(NotPeriodicError):
            test_cell = topology_graph.get_periodic_cell()

    else:
        test_cell = topology_graph.get_periodic_cell()
        assert np.all(np.array([
            np.allclose(i, j, atol=1e-4)
            for i, j in zip(test_cell, periodic_cell)
        ]))
