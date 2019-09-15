import pytest
import stk
import numpy as np


@pytest.fixture(scope='function')
def tmp_vertex():
    vdata = stk.VertexData(1, 2, 3)
    vdata.cell = np.array([3, 4, 12])
    vdata.id = 10
    e1, e2, e3 = (
        stk.EdgeData(vdata), stk.EdgeData(vdata), stk.EdgeData(vdata)
    )
    e1.id = 12
    e2.id = 24
    e3.id = 36
    return stk.Vertex(vdata)


@pytest.fixture(scope='session')
def vertex():
    vdata = stk.VertexData(1, 2, 3)
    vdata.cell = np.array([3, 4, 12])
    vdata.id = 10
    e1, e2, e3 = (
        stk.EdgeData(vdata), stk.EdgeData(vdata), stk.EdgeData(vdata)
    )
    e1.id = 12
    e2.id = 24
    e3.id = 36
    return stk.Vertex(vdata)


@pytest.fixture(scope='session')
def graph_components():
    vdata1, vdata2, vdata3 = vertex_data = (
        stk.VertexData(0, 0, 0),
        stk.VertexData(-1, -2, -3),
        stk.VertexData(1, 2, 3)
    )
    for i, vertex in enumerate(vertex_data):
        vertex.id = i
    vertices = tuple(stk.Vertex(data) for data in vertex_data)

    edges = (
        stk.EdgeData(vdata1, vdata2).get_edge(),
        stk.EdgeData(vdata1, vdata3).get_edge()
    )
    return vertices, edges


@pytest.fixture(scope='session')
def graph_components_alt1():
    vdata1, vdata2, vdata3, vdata4 = vertex_data = (
        stk.VertexData(0, 0, 0),
        stk.VertexData(1, 0, 0),
        stk.VertexData(-1, 0, 0),
        stk.VertexData(0, 1, 0)
    )
    for i, vertex in enumerate(vertex_data):
        vertex.id = i
    vertices = tuple(stk.Vertex(data) for data in vertex_data)

    edges = (
        stk.EdgeData(vdata1, vdata2).get_edge(),
        stk.EdgeData(vdata1, vdata3).get_edge(),
        stk.EdgeData(vdata1, vdata4).get_edge()
    )
    return vertices, edges


@pytest.fixture(scope='session')
def periodic_graph_components():
    vdata1, vdata2, vdata3, vdata4 = vertex_data = (
        stk.VertexData(0, 0, 0),
        stk.VertexData(1, 0, 0),
        stk.VertexData(-1, 0, 0),
        stk.VertexData(0, 1, 0)
    )
    for i, vertex in enumerate(vertex_data):
        vertex.id = i
    vertices = tuple(stk.Vertex(data) for data in vertex_data)

    edges = (
        stk.EdgeData(vdata1, vdata2).get_edge(),
        stk.EdgeData(
            vdata1,
            vdata3,
            periodicity=[1, 0, 0],
            lattice_constants=(
                np.array([1, 0, 0]),
                np.array([0, 1, 0]),
                np.array([0, 0, 1])
            )
        ).get_edge()
    )
    return vertices, edges


@pytest.fixture(scope='session')
def periodic_graph_components_alt1():
    vdata1, vdata2, vdata3, vdata4 = vertex_data = (
        stk.VertexData(0, 0, 0),
        stk.VertexData(1, 0, 0),
        stk.VertexData(-1, 0, 0),
        stk.VertexData(0, 1, 0)
    )
    for i, vertex in enumerate(vertex_data):
        vertex.id = i
    vertices = tuple(stk.Vertex(data) for data in vertex_data)

    edges = (
        stk.EdgeData(vdata1, vdata2).get_edge(),
        stk.EdgeData(vdata1, vdata3).get_edge(),
        stk.EdgeData(
            vdata1,
            vdata4,
            periodicity=[0, 0, 1],
            lattice_constants=(
                np.array([1, 0, 0]),
                np.array([0, 1, 0]),
                np.array([0, 0, 1])
            )
        ).get_edge()
    )
    return vertices, edges


@pytest.fixture(scope='session')
def ab_chain3():
    return stk.polymer.Linear('AB', 3)


@pytest.fixture(scope='function')
def tmp_ab_chain3():
    return stk.polymer.Linear('AB', 3)


@pytest.fixture(scope='session')
def ab_chain6():
    return stk.polymer.Linear('AB', 6)


@pytest.fixture(scope='session')
def honeycomb_lattice():
    return stk.cof.Honeycomb((3, 3, 1))
