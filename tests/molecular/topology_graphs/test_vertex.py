import stk


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
