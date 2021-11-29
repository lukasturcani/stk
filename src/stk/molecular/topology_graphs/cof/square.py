"""
Square
======

"""

import numpy as np

from ..topology_graph import Edge
from .cof import Cof
from .vertices import LinearVertex, NonLinearVertex


class Square(Cof):
    """
    Represents a sqaure COF topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.PeriodicSquare(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCC(Br)',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)C(Br)=C1Br',
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
                if all(p == 0 for p in bond.get_periodicity())
            ),
        )

    ``Collapser(scale_steps=False)`` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.PeriodicSquare(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCC(Br)',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrC1=C(Br)C(Br)=C1Br',
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
                if all(p == 0 for p in bond.get_periodicity())
            ),
        )

    Building blocks with four and two functional groups are required
    for this topology graph.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cof-topology-graph-examples`:
    *Multi-Building Block COF Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0
        | 2-functional groups: 1 to 2

    See :class:`.Cof` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0., 1., 0.]),
        np.array([0., 0., 1.])
    )

    _non_linears = (
        NonLinearVertex(0, (0.5)*_a + (0.5)*_b + (0.5)*_c),
    )
    _vertex_prototypes = (
        *_non_linears,
        LinearVertex.init_at_shifted_center(
            id=1,
            vertices=(_non_linears[0], _non_linears[0]),
            cell_shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants,
        ),
        LinearVertex.init_at_shifted_center(
            id=2,
            vertices=(_non_linears[0], _non_linears[0]),
            cell_shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants,
        ),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[1], _vertex_prototypes[0]),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[1],
            vertex2=_vertex_prototypes[0],
            periodicity=(1, 0, 0),
        ),
        Edge(2, _vertex_prototypes[2], _vertex_prototypes[0]),
        Edge(
            id=3,
            vertex1=_vertex_prototypes[2],
            vertex2=_vertex_prototypes[0],
            periodicity=(0, 1, 0),
        ),
    )
