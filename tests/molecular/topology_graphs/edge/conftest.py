import pytest
import stk


class CaseData:
    def __init__(
        self,
        edge,
        id,
        vertex1_id,
        vertex2_id,
        periodicity,
        is_periodic,
    ):
        self.edge = edge
        self.id = id
        self.vertex1_id = vertex1_id
        self.vertex2_id = vertex2_id
        self.periodicity = periodicity
        self.is_periodic = is_periodic


@pytest.fixture
def case_data(id, vertex1, vertex2, periodicity):
    return CaseData(
        edge=stk.Edge(
            id=id,
            vertex1=vertex1,
            vertex2=vertex2,
            periodicity=periodicity,
        ),
        id=id,
        vertex1_id=vertex1.get_id(),
        vertex2_id=vertex2.get_id(),
        periodicity=periodicity,
        is_periodic=any(i != 0 for i in periodicity),
    )


@pytest.fixture(
    params=(
        stk.Vertex(0, [1, 2, 3]),
        stk.Vertex(3, [4, 2, 1]),
    ),
)
def vertex(request):
    return request.param


@pytest.fixture
def vertex1(vertex):
    return vertex.clone()


@pytest.fixture
def vertex2(vertex):
    return vertex.clone()


@pytest.fixture(
    params=(
        (0, 0, 0),
        (0, 1, 0),
        (-1, 1, 1),
        (20, -1, 12),
    ),
)
def periodicity(request):
    return request.param


@pytest.fixture(params=(0, 12))
def id(request):
    return request.param
