import numpy as np


def test_placement(case_data):
    """
    Test :meth:`.Vertex.place_building_block`.

    The point of this test is to place a building block on two
    vertices, one which flips and one which does not, to guarantee
    that they produce the opposite result.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the vertices to test and the building block
        to place on them.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_placement(
        vertex1=case_data.vertex1,
        vertex2=case_data.vertex2,
        building_block=case_data.building_block,
        atom_ids=case_data.atom_ids,
    )


def _test_placement(vertex1, vertex2, building_block, atom_ids):
    """
    Test :meth:`.Vertex.place_building_block`.

    The point of this test is to place a building block on two
    vertices, one which flips and one which does not, to guarantee
    that they produce the opposite result.

    Parameters
    ----------
    vertex1 : :class:`.Vertex`
        The first vertex to test.

    vertex2 : :class:`.Vertex`
        The second vertex to test.

    building_block : :class:`.BuildingBlock`
        The building block to place.

    atom_ids : :class:`tuple` of :class:`int`
        The ids of atoms used for calculating the building block
        orientation.

    Returns
    -------
    None : :class:`NoneType`


    """

    position_matrix1 = vertex1.place_building_block(building_block, ())
    building_block1 = building_block.with_position_matrix(
        position_matrix=position_matrix1,
    )
    position_matrix2 = vertex2.place_building_block(building_block, ())
    building_block2 = building_block.with_position_matrix(
        position_matrix=position_matrix2,
    )
    normal1 = building_block1.get_plane_normal(
        atom_ids=atom_ids,
    )
    normal2 = building_block2.get_plane_normal(
        atom_ids,
    )
    assert np.allclose(
        a=normal1 @ normal2,
        b=[-1, 0, 0],
        atol=1e-13,
    )
