import stk
import numpy as np


def test_init_vertex_data():
    data = stk.VertexData(1, 2, 3)
    assert all(data.position == [1, 2, 3])
    assert data.edges == []
    assert all(data.cell == [0, 0, 0])


def test_init_at_center():
    vertices = [stk.VertexData(1, 2, 3) for i in range(10)]
    data = stk.VertexData.init_at_center(*vertices)
    assert all(data.position == [1, 2, 3])
    assert data.edges == []
    assert all(data.cell == [0, 0, 0])

    vertices = (stk.VertexData(-1, -2, -3), stk.VertexData(1, 2, 3))
    data = stk.VertexData.init_at_center(*vertices)
    assert all(data.position == [0, 0, 0])
    assert data.edges == []
    assert all(data.cell == [0, 0, 0])


def test_init_vertex():
    data = stk.VertexData(1, 2, 3, np.array([10, 4, 2]))
    vertex = stk.Vertex(12, data)

    assert vertex.id == 12
    assert all(vertex.get_position() == [1, 2, 3])
    assert list(vertex.get_edge_ids()) == [12, 24, 36]
    assert all(vertex.get_cell() == [10, 4, 2])

    data = stk.VertexData(10, 20, 30)
    vertex = stk.Vertex(21, data)

    assert vertex.id == 21
    assert all(vertex.get_position() == [10, 20, 30])
    assert list(vertex.get_edge_ids()) == []
    assert all(vertex.get_cell() == [0, 0, 0])
