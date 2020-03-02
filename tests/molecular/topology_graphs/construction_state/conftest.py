import pytest
import stk


class _TestCase:
    def __init__(
        self,
        construction_state,
        atom_infos,
        atoms,
        bond_infos,
        bonds,
        building_blocks,
        building_block_counts,
        edges,
        edge_group,
        functional_groups,
        vertex_edges,
        num_edges,
        num_vertices,
        position_matrix,
        vertices,
    ):

        self.construction_state = construction_state
        self.atom_infos = atom_infos
        self.atoms = atoms
        self.bond_infos = bond_infos
        self.bonds = bonds
        self.building_blocks = building_blocks
        self.building_block_counts = building_block_counts
        self.edges = edges
        self.edge_group = edge_group
        self.functional_group = functional_groups
        self.vertex_edges = vertex_edges
        self.num_edges = num_edges
        self.num_vertices = num_vertices
        self.position_matrix = position_matrix
        self.vertices = vertices


@pytest.fixture
def test_case(building_block_vertices):
    return _TestCase(
        construction_state=stk.ConstructionState(
            building_block_vertices=building_block_vertices,
            edges=get_edges(building_block_vertices),
            scale=1,
        ),
    )


def get_edges(building_block_vertices):
    vertices = (
        vertex
        for vertices in building_block_vertices.values()
        for vertex in vertices
    )
    vertex1 = next(vertices)
    for id_, vertex2 in enumerate(vertices):
        yield stk.Edge(id_, vertex1, vertex2)
        vertex1 = vertex2


@pytest.fixture
def building_block_vertices(building_block1, building_block2):
    return {
        building_block1: (
            stk.Vertex(0, [0, 0, 0]),
            stk.Vertex(1, [10, 0, 0]),
        ),
        building_block2: (
            stk.Vertex(2, [20, 0, 0]),
            stk.Vertex(3, [30, 0, 0]),
        ),
    }


@pytest.fixture(
    params=(
        stk.BuildingBlock('BrCC', [stk.BromoFactory()]),
    ),
)
def building_block1(request):
    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
    ),
)
def building_block2(request):
    return request.param
