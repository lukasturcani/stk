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
    vdata = stk.VertexData(1, 2, 3, np.array([10, 4, 2]))
    e1, e2, e3 = (
        stk.EdgeData(vdata), stk.EdgeData(vdata), stk.EdgeData(vdata)
    )
    e1.id = 12
    e2.id = 24
    e3.id = 36
    vertex = stk.Vertex(vdata)

    assert vertex.id is None
    assert all(vertex.get_position() == [1, 2, 3])
    assert list(vertex.get_edge_ids()) == [12, 24, 36]
    assert all(vertex.get_cell() == [10, 4, 2])

    vdata = stk.VertexData(10, 20, 30)
    e1, e2 = stk.EdgeData(vdata), stk.EdgeData(vdata)
    e1.id = 100
    e2.id = 200
    vertex = stk.Vertex(vdata)

    assert vertex.id is None
    assert all(vertex.get_position() == [10, 20, 30])
    assert list(vertex.get_edge_ids()) == [100, 200]
    assert all(vertex.get_cell() == [0, 0, 0])
