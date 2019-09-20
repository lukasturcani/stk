import stk


def test_reference_creation():
    v1, v2, v3 = vertices = (
        stk.VertexData(1, 2, 3),
        stk.VertexData(10, 20, 30),
        stk.VertexData(100, 200, 300)
    )

    e1 = stk.EdgeData(v1, v2)
    e2 = stk.EdgeData(v1, v2, v3)

    for vertex in vertices[:2]:
        assert vertex in e1.vertices
        assert e1 in vertex.edges
        assert len(vertex.edges) == 2
    assert len(v3.edges) == 1

    for vertex in vertices:
        assert vertex in e2.vertices
        assert e2 in vertex.edges
