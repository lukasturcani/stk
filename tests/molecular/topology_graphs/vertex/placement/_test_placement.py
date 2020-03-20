import numpy as np


def test_placement(case_data):
    """
    Test placement methods.

    This tests both :meth:`.Vertex.place_building_block` and
    :meth:`.Vertex.map_functional_groups_to_edges`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the vertex to test, the building block which
        it  places, and the correct position and alignment of the
        building block after placement.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_placement(
        # Use a clone to ensure that the clone() method is implemented
        # correctly.
        vertex=case_data.vertex.clone(),
        edges=case_data.edges,
        building_block=case_data.building_block,
        position=case_data.position,
        alignment_tests=case_data.alignment_tests,
        functional_group_edges=case_data.functional_group_edges,
        position_ids=case_data.position_ids,
    )


def _test_placement(
    vertex,
    edges,
    building_block,
    position,
    alignment_tests,
    functional_group_edges,
    position_ids,
):
    """
    Test placement methods.

    This tests both :meth:`.Vertex.place_building_block` and
    :meth:`.Vertex.map_functional_groups_to_edges`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    edges : :class:`tuple` of :class:`.Edge`
        The edges connected to `vertex`.

    building_block : :class:`.BuildingBlock`
        The building block which is placed by `vertex`.

    position : :class:`numpy.ndarray`
        The correct centroid of atoms in `position_ids` after
        placement.

    alignment_tests : :class:`dict`
        Maps a :class:`callable` to a :class:`numpy.ndarray`.
        Each :class:`callable` takes a singe parameter,
        `building_block`, and returns a :class:`numpy.ndarray`.

        Each :class:`callable` in `alignment_tests` represents a
        test, and the array which it maps to is the correct result
        for that test. If the array returned by the :class:`callable`
        does not match the array it is mapped to, the test will fail.

        The point of these tests is to make sure that `building_block`
        is aligned correctly after placement. Therefore, alignment
        tests should return some vector which depends on the specific
        orientation of the building block.

    functional_group_edges : :class:`dict`
        The correct mapping of functional group id to edge id.

    position_ids : :class:`tuple` of :class:`int`
        The ids of atoms which are used to determine if the building
        block was positioned correctly.

    Returns
    -------
    None : :class:`NoneType`

    """

    position_matrix = vertex.place_building_block(
        building_block=building_block,
        edges=edges,
    )
    building_block = building_block.with_position_matrix(
        position_matrix=position_matrix,
    )
    assert np.allclose(
        a=building_block.get_centroid(position_ids),
        b=position,
        atol=1e-14,
    )
    for test, result in alignment_tests.items():
        # Do the assert via a function call here so that pytest prints
        # the values being compared when an error occurs.
        assert_equal_vectors(test(building_block), result)

    functional_group_edges_ = vertex.map_functional_groups_to_edges(
        building_block=building_block,
        edges=edges,
    )
    assert functional_group_edges_ == functional_group_edges


def assert_equal_vectors(vector1, vector2):
    assert np.all(np.equal(vector2, vector1))
