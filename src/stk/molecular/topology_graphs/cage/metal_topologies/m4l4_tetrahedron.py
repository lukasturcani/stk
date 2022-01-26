"""
M4L4 Tetrahedron
================

"""

import numpy as np

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import NonLinearVertex


class M4L4Tetrahedron(Cage):
    """
    Represents a cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        iron_atom = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
            smiles='C1=NC(C=NBr)=CC=C1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralDelta(
                metals=iron_atom,
                ligands=bb2,
                optimizer=stk.MCHammer(),
            ),
        )

        # Assign Bromo functional groups to the metal complex.
        iron_oct_delta = stk.BuildingBlock.init_from_molecule(
            molecule=complex,
            functional_groups=[stk.BromoFactory()],
        )

        # Define building blocks.
        bb3 = stk.BuildingBlock(
            smiles=(
                'C1=C(C=C(C=C1Br)Br)Br'
            ),
            functional_groups=[stk.BromoFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L4Tetrahedron(
                building_blocks={
                    iron_oct_delta: (0, 1, 2, 3),
                    bb3: (4, 5, 6, 7),
                },
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cage.get_atoms(),
                    cage.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        1
                        if bond.get_order() == 9
                        else bond.get_order()
                    ),
                ) for bond in cage.get_bonds()
            ),
        )

    :class:`.MCHammer` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        iron_atom = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
            smiles='C1=NC(C=NBr)=CC=C1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralDelta(
                metals=iron_atom,
                ligands=bb2,
                optimizer=stk.MCHammer(),
            ),
        )

        # Assign Bromo functional groups to the metal complex.
        iron_oct_delta = stk.BuildingBlock.init_from_molecule(
            molecule=complex,
            functional_groups=[stk.BromoFactory()],
        )

        # Define building blocks.
        bb3 = stk.BuildingBlock(
            smiles=(
                'C1=C(C=C(C=C1Br)Br)Br'
            ),
            functional_groups=[stk.BromoFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L4Tetrahedron(
                building_blocks={
                    iron_oct_delta: (0, 1, 2, 3),
                    bb3: (4, 5, 6, 7),
                },
                optimizer=stk.MCHammer(),
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cage.get_atoms(),
                    cage.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        1
                        if bond.get_order() == 9
                        else bond.get_order()
                    ),
                ) for bond in cage.get_bonds()
            ),
        )

    Building blocks with three functional groups are required for this
    topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups (metal): 0 to 3
        | 3-functional groups (linker): 4 to 7

    See :class:`.Cage` for more details and examples.

    """

    _non_linears = (
        NonLinearVertex(0, [0, 0, np.sqrt(6)/2]),
        NonLinearVertex(1, [-1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        NonLinearVertex(2, [1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        NonLinearVertex(3, [0, 2*np.sqrt(3)/3, -np.sqrt(6)/6]),
    )

    _vertex_prototypes = (
        *_non_linears,

        NonLinearVertex.init_at_center(
            id=4,
            vertices=(
                _non_linears[0],
                _non_linears[1],
                _non_linears[2],
            ),
        ),
        NonLinearVertex.init_at_center(
            id=5,
            vertices=(
                _non_linears[0],
                _non_linears[1],
                _non_linears[3],
            ),
        ),
        NonLinearVertex.init_at_center(
            id=6,
            vertices=(
                _non_linears[0],
                _non_linears[2],
                _non_linears[3],
            ),
        ),
        NonLinearVertex.init_at_center(
            id=7,
            vertices=(
                _non_linears[1],
                _non_linears[2],
                _non_linears[3],
            ),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[7]),
        Edge(6, _vertex_prototypes[2], _vertex_prototypes[4]),
        Edge(7, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(8, _vertex_prototypes[2], _vertex_prototypes[7]),
        Edge(9, _vertex_prototypes[3], _vertex_prototypes[5]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[6]),
        Edge(11, _vertex_prototypes[3], _vertex_prototypes[7]),
    )

    _num_windows = 4
    _num_window_types = 1
