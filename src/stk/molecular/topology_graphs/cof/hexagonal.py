"""
Hexagonal
=========

"""

import numpy as np

from ..topology_graph import Edge
from .cof import Cof
from .vertices import LinearVertex, NonLinearVertex


class Hexagonal(Cof):
    """
    Represents a hexagonal COF topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='Brc1c(Br)c(Br)c(Br)c(Br)c1Br',
            functional_groups=[stk.BromoFactory()],
        )

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.Hexagonal(
                building_blocks=(bb1, bb2),
                lattice_size=(2, 2, 1),
            ),
        )
        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cof.get_atoms(),
                    cof.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in cof.get_bonds()
            ),
        )

    ``Collapser(scale_steps=False)`` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='Brc1c(Br)c(Br)c(Br)c(Br)c1Br',
            functional_groups=[stk.BromoFactory()],
        )

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.Hexagonal(
                building_blocks=(bb1, bb2),
                lattice_size=(2, 2, 1),
                optimizer=stk.Collapser(scale_steps=False),
            ),
        )
        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cof.get_atoms(),
                    cof.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in cof.get_bonds()
            ),
        )

    Building blocks with six and two functional groups are required
    for this topology graph.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cof-topology-graph-examples`:
    *Multi-Building Block COF Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 6-functional groups: 0 to 3
        | 2-functional groups: 4 to 15

    See :class:`.Cof` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0]),
        np.array([0, 0, 5/1.7321])
    )

    _non_linears = (
        NonLinearVertex(0, (1/4)*_a + (1/4)*_b + (1/2)*_c),
        NonLinearVertex(1, (1/4)*_a + (3/4)*_b + (1/2)*_c),
        NonLinearVertex(2, (3/4)*_a + (1/4)*_b + (1/2)*_c),
        NonLinearVertex(3, (3/4)*_a + (3/4)*_b + (1/2)*_c),
    )

    _vertex_prototypes = (
        *_non_linears,
        LinearVertex.init_at_center(
            id=4,
            vertices=(_non_linears[0], _non_linears[1]),
        ),
        LinearVertex.init_at_center(
            id=5,
            vertices=(_non_linears[0], _non_linears[2]),
        ),
        LinearVertex.init_at_center(
            id=6,
            vertices=(_non_linears[1], _non_linears[2]),
        ),
        LinearVertex.init_at_center(
            id=7,
            vertices=(_non_linears[1], _non_linears[3]),
        ),
        LinearVertex.init_at_center(
            id=8,
            vertices=(_non_linears[2], _non_linears[3]),
        ),
        LinearVertex.init_at_shifted_center(
            id=9,
            vertices=(_non_linears[0], _non_linears[2]),
            cell_shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants,
        ),
        LinearVertex.init_at_shifted_center(
            id=10,
            vertices=(_non_linears[0], _non_linears[1]),
            cell_shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        LinearVertex.init_at_shifted_center(
            id=11,
            vertices=(_non_linears[0], _non_linears[3]),
            cell_shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        LinearVertex.init_at_shifted_center(
            id=12,
            vertices=(_non_linears[2], _non_linears[1]),
            cell_shifts=((0, 0, 0), (1, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        LinearVertex.init_at_shifted_center(
            id=13,
            vertices=(_non_linears[2], _non_linears[3]),
            cell_shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        LinearVertex.init_at_shifted_center(
            id=14,
            vertices=(_non_linears[1], _non_linears[3]),
            cell_shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants,
        ),
        LinearVertex.init_at_shifted_center(
            id=15,
            vertices=(_non_linears[3], _non_linears[0]),
            cell_shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants,
        )
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[4], _vertex_prototypes[1]),

        Edge(2, _vertex_prototypes[5], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[5], _vertex_prototypes[2]),

        Edge(4, _vertex_prototypes[6], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[6], _vertex_prototypes[2]),

        Edge(6, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[7], _vertex_prototypes[3]),

        Edge(8, _vertex_prototypes[8], _vertex_prototypes[2]),
        Edge(9, _vertex_prototypes[8], _vertex_prototypes[3]),

        Edge(10, _vertex_prototypes[9], _vertex_prototypes[0]),
        Edge(
            id=11,
            vertex1=_vertex_prototypes[9],
            vertex2=_vertex_prototypes[2],
            periodicity=(-1, 0, 0),
        ),

        Edge(12, _vertex_prototypes[10], _vertex_prototypes[0]),
        Edge(
            id=13,
            vertex1=_vertex_prototypes[10],
            vertex2=_vertex_prototypes[1],
            periodicity=(0, -1, 0),
        ),

        Edge(14, _vertex_prototypes[11], _vertex_prototypes[0]),
        Edge(
            id=15,
            vertex1=_vertex_prototypes[11],
            vertex2=_vertex_prototypes[3],
            periodicity=(0, -1, 0),
        ),

        Edge(16, _vertex_prototypes[12], _vertex_prototypes[2]),
        Edge(
            id=17,
            vertex1=_vertex_prototypes[12],
            vertex2=_vertex_prototypes[1],
            periodicity=(1, -1, 0),
        ),

        Edge(18, _vertex_prototypes[13], _vertex_prototypes[2]),
        Edge(
            id=19,
            vertex1=_vertex_prototypes[13],
            vertex2=_vertex_prototypes[3],
            periodicity=(0, -1, 0),
        ),

        Edge(20, _vertex_prototypes[14], _vertex_prototypes[1]),
        Edge(
            id=21,
            vertex1=_vertex_prototypes[14],
            vertex2=_vertex_prototypes[3],
            periodicity=(-1, 0, 0),
        ),

        Edge(22, _vertex_prototypes[15], _vertex_prototypes[3]),
        Edge(
            id=23,
            vertex1=_vertex_prototypes[15],
            vertex2=_vertex_prototypes[0],
            periodicity=(1, 0, 0),
        ),
    )
