"""
Kagome
======

"""

import numpy as np

from ..topology_graph import Edge
from .cof import Cof
from .vertices import LinearVertex, NonLinearVertex


class Kagome(Cof):
    """
    Represents a kagome COF topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.Kagome(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCC(Br)',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrC1C(Br)CC(Br)C(Br)C1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                lattice_size=(3, 3, 1),
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

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.Kagome(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCC(Br)',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrC1C(Br)CC(Br)C(Br)C1',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                lattice_size=(3, 3, 1),
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


    Building blocks with four and two functional groups are required
    for this topology graph.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cof-topology-graph-examples`:
    *Multi-Building Block COF Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0 to 2
        | 2-functional groups: 3 to 8

    See :class:`.Cof` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0.]),
        np.array([0., 0., 5/1.7321])
    )

    _non_linears = (
        NonLinearVertex(0, (1/4)*_a + (3/4)*_b + (0.5)*_c),
        NonLinearVertex(1, (3/4)*_a + (3/4)*_b + (1/2)*_c),
        NonLinearVertex(2, (3/4)*_a + (1/4)*_b + (1/2)*_c),
    )

    _vertex_prototypes = (
        *_non_linears,
        LinearVertex.init_at_center(
            id=3,
            vertices=(_non_linears[0], _non_linears[1]),
        ),
        LinearVertex.init_at_center(
            id=4,
            vertices=(_non_linears[0], _non_linears[2]),
        ),
        LinearVertex.init_at_center(
            id=5,
            vertices=(_non_linears[1], _non_linears[2]),
        ),
        LinearVertex.init_at_shifted_center(
            id=6,
            vertices=(_non_linears[0], _non_linears[1]),
            cell_shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        LinearVertex.init_at_shifted_center(
            id=7,
            vertices=(_non_linears[0], _non_linears[2]),
            cell_shifts=((0, 0, 0), (-1, 1, 0)),
            lattice_constants=_lattice_constants
        ),
        LinearVertex.init_at_shifted_center(
            id=8,
            vertices=(_non_linears[1], _non_linears[2]),
            cell_shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants
        ),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[3], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[3], _vertex_prototypes[1]),

        Edge(2, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[4], _vertex_prototypes[2]),

        Edge(4, _vertex_prototypes[5], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[5], _vertex_prototypes[2]),

        Edge(6, _vertex_prototypes[6], _vertex_prototypes[0]),
        Edge(
            id=7,
            vertex1=_vertex_prototypes[6],
            vertex2=_vertex_prototypes[1],
            periodicity=(-1, 0, 0),
        ),

        Edge(8, _vertex_prototypes[7], _vertex_prototypes[0]),
        Edge(
            id=9,
            vertex1=_vertex_prototypes[7],
            vertex2=_vertex_prototypes[2],
            periodicity=(-1, 1, 0),
        ),

        Edge(10, _vertex_prototypes[8], _vertex_prototypes[1]),
        Edge(
            id=11,
            vertex1=_vertex_prototypes[8],
            vertex2=_vertex_prototypes[2],
            periodicity=(0, 1, 0),
        ),
    )
