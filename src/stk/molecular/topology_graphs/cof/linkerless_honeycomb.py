"""
Linkerless Honeycomb
====================

"""

import numpy as np

from ..topology_graph import Edge
from .cof import Cof
from .vertices import NonLinearVertex


class LinkerlessHoneycomb(Cof):
    """
    Represents a honeycomb COF topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.LinkerlessHoneycomb(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCC(CBr)CBr',
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
            topology_graph=stk.cof.LinkerlessHoneycomb(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCC(CBr)CBr',
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


    Building blocks with three functional groups are required
    for this topology graph.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cof-topology-graph-examples`:
    *Multi-Building Block COF Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 1

    See :class:`.Cof` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0.]),
        np.array([0., 0., 5/1.7321]),
    )

    _vertex_prototypes = (
        NonLinearVertex(0, (1/3)*_a + (1/3)*_b + (1/2)*_c),
        NonLinearVertex(1, (2/3)*_a + (2/3)*_b + (1/2)*_c),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[1]),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            periodicity=(-1, 0, 0),
        ),
        Edge(
            id=2,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            periodicity=(0, -1, 0),
        ),
    )
