import numpy as np


def test_apply_scale(make_vertex, position, scale):
    vertex = make_vertex(position)
    before = vertex.get_position()
    vertex.apply_scale(scale)
    assert np.allclose(
        a=vertex.get_position(),
        b=scale*before,
        atol=1e-32,
    )


def test_clone(make_vertex):
    vertex = make_vertex()
    vertex.attr = 1
    vertex._attr = 2
    clone = vertex.clone()
    assert clone.get_position() == vertex.get_position()
    assert clone.attr == vertex.attr
    assert not hasattr(clone, '_attr')
    assert clone.get_num_edges() == vertex.get_num_edges()
    assert list(clone.get_edge_ids()) == list(vertex.get_edge_ids())
    assert clone.get_cell() == vertex.get_cell()


def test_get_position(make_vertex, position):
    vertex = make_vertex(position)
    assert np.all(np.equal(vertex.get_position(), position))


def test_get_num_edges(make_vertex, make_edges, edge_ids):
    vertex = make_vertex(edges=make_edges(edge_ids))
    assert vertex.get_num_edges() == len(edge_ids)


def test_get_edge_ids(make_vertex, make_edges, edge_ids):
    vertex = make_vertex(edges=make_edges(edge_ids))
    assert tuple(vertex.get_edge_ids()) == edge_ids


def test_get_cell(make_vertex, cell):
    vertex = make_vertex(cell=cell)
    assert vertex.get_cell() == cell


def test_set_constructed_molecule(make_vertex, constructed_molecule):
    vertex = make_vertex()
    vertex.set_constructed_molecule(molecule)
    assert


class TestPlacement:
    def _test_placement(self, test_case, building_block, position):
        test_case.vertex.place_building_block(
            building_block=building_block,
            vertices=test_case.vertices,
            edges=test_case.edges,
        )
        assert (
            np.all(np.equal(building_block.get_centroid(), position))
        )

    def _test_alignment(self, test_case, building_block):
        bonder_centroids = enumerate(
            building_block.get_bonder_centroids()
        )
        for i, bonder_centroid in bonder_centroids:
            closest_edge = min(
                test_case.vertex.get_edge_ids(),
                key=lambda edge_id: distance(
                    position1=bonder_centroid,
                    position2=test_case.edges[edge_id].get_position(),
                ),
            )
            assert closest_edge == test_case.alignments[i]

    def _test_assignment(self, test_case, building_block):
        vertex = test_case.vertex
        functional_group_edges = vertex.assign_func_groups_to_edges(
            building_block=building_block,
            vertices=test_case.vertices,
            edges=test_case.edges,
        )
        assert (
            functional_group_edges == test_case.functional_group_edges
        )

    def _test_place_building_block(
        self,
        make_test_case,
        building_block,
        position,
    ):
        test_case = make_test_case(position)
        self._test_placement(test_case, building_block, position)
        self._test_alignment(test_case, building_block)
        self._test_assignment(test_case, building_block)

    def test_place_building_block_0(
        self,
        make_test_case_0,
        building_block_0,
        position,
    ):
        _test_place_building_block(
            make_test_case=make_test_case_0,
            building_block=building_block_0,
            position=position,
        )

    def test_place_building_block_1(
        self,
        make_test_case_1,
        building_block_1,
        position,
    ):
        _test_place_building_block(
            make_test_case=make_test_case_1,
            building_block=building_block_1,
            position=position,
        )

    def test_place_building_block_2(
        self,
        make_test_case_2,
        building_block_2,
        position,
    ):
        _test_place_building_block(
            make_test_case=make_test_case_2,
            building_block=building_block_2,
            position=position,
        )

    def test_place_building_block_3(
        self,
        make_test_case_3,
        building_block_3,
        position,
    ):
        _test_place_building_block(
            make_test_case=make_test_case_3,
            building_block=building_block_3,
            position=position,
        )

    def test_place_building_block_4(
        self,
        make_test_case_4,
        building_block_4,
        position,
    ):
        _test_place_building_block(
            make_test_case=make_test_case_4,
            building_block=building_block_4,
            position=position,
        )

    def test_place_building_block_5(
        self,
        make_test_case_5,
        building_block_5,
        position,
    ):

        _test_place_building_block(
            make_test_case=make_test_case_5,
            building_block=building_block_5,
            position=position,
        )


def test_get_edge_centroid(vertex_test_case, edges):
    vertex = vertex_test_case.vertex
    centroid_edges = vertex_test_case.centroid_edges
    vertices = vertex_test_case.vertices
    assert (
        vertex._get_edge_centroid(centroid_edges, vertices)
        == vertex_test_case.edge_centroid
    )


def test_get_edge_plane_normal(vertex_test_case, edges):
    vertex = vertex_test_case.vertex


def test_get_molecule_centroid(vertex_test_case, edges, molecule):
    assert False
