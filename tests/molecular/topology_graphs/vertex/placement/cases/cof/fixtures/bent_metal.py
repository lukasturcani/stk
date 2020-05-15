import pytest
import stk
from functools import partial
from scipy.spatial.distance import euclidean

from ....case_data import CaseData

vertices = stk.cage.vertices

metal_atom = stk.BuildingBlock(
    smiles='[Pd+2]',
    functional_groups=(
        stk.SingleAtom(stk.Pd(0, charge=2))
        for i in range(4)
    ),
    position_matrix=([0, 0, 0], ),
)

palladium_bi_1 = stk.BuildingBlock(
    smiles='[H]N([H])C([H])([H])C([H])([H])N([H])[H]',
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#7]~[#6]',
            bonders=(0, ),
            deleters=(),
        ),
    ]
)
palladium_cispbi_sqpl = stk.ConstructedMolecule(
    stk.metal_complex.CisProtectedSquarePlanar(
        metals={metal_atom: 0},
        ligands={palladium_bi_1: 0},
        reaction_factory=stk.DativeReactionFactory(
            stk.GenericReactionFactory(
                bond_orders={
                    frozenset({
                        stk.GenericFunctionalGroup,
                        stk.SingleAtom
                    }): 9
                }
            )
        )
    )
)


@pytest.fixture
def bent_metal(position, bent_aligner_edge, bent_building_block):

    point1, point2 = points = (
        position + [10, 0, 0],
        position + [-10, 0, 0],
    )

    vertex = vertices._BentMetalComplexCageVertex(
        id=0,
        position=position,
        aligner_edge=bent_aligner_edge,
    )
    edges = tuple(get_bent_edges(vertex))

    def get_point_closest_to_edge_centroid(building_block):
        return get_closest_point(
            points=points,
            point=get_edge_centroid_position(edges),
        )

    def get_point_closest_to_core(building_block):
        return get_closest_point(
            points=points,
            point=get_core_position(building_block),
        )

    return CaseData(
        vertex=vertex,
        edges=edges,
        building_block=bent_building_block,
        position=position,
        alignment_tests={
            get_point_closest_to_edge_centroid: point1,
            get_point_closest_to_core: point2,
        },
        functional_group_edges=({0: 0, 1: 1}),
    )


def get_edge_centroid_position(edges):
    edge_centroid = (
        sum(edge.get_position() for edge in edges) / len(edges)
    )
    return edge_centroid


def get_core_position(building_block):
    return building_block.get_centroid(
        atom_ids=building_block.get_core_atom_ids(),
    )


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


def get_bent_edges(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [10, 10, 0])
    vertex3 = stk.Vertex(2, vertex.get_position() + [10, -10, 0])
    yield stk.Edge(0, vertex, vertex2)
    yield stk.Edge(1, vertex, vertex3)


@pytest.fixture(params=(0, 1))
def bent_aligner_edge(request):
    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock.init_from_molecule(
            molecule=palladium_cispbi_sqpl,
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[Pd]~[#7]',
                    bonders=(0, ),
                    deleters=(),
                ),
            ]
        ),
    ),
)
def bent_building_block(request):
    return request.param
