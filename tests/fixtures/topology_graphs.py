import pytest
import stk
import numpy as np


@pytest.fixture(scope='function')
def tmp_vertex():
    vdata = stk.VertexData(1, 2, 3, np.array(3, 4, 12))
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
    vdata = stk.VertexData(1, 2, 3, np.array(3, 4, 12))
    vdata.id = 10
    e1, e2, e3 = (
        stk.EdgeData(vdata), stk.EdgeData(vdata), stk.EdgeData(vdata)
    )
    e1.id = 12
    e2.id = 24
    e3.id = 36
    return stk.Vertex(vdata)


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
