"""
M4L6 Tetrahedron with Spacer
============================

"""

import numpy as np

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex


class M4L6TetrahedronSpacer(Cage):
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
                'C1=CC(=CC=C1C2=CC=C(C=C2)Br)Br'
            ),
            functional_groups=[stk.BromoFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L6TetrahedronSpacer(
                building_blocks={
                    iron_oct_delta: (0, 1, 2, 3),
                    bb3: (4, 5, 6, 7, 8, 9),
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
                'C1=CC(=CC=C1C2=CC=C(C=C2)Br)Br'
            ),
            functional_groups=[stk.BromoFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L6TetrahedronSpacer(
                building_blocks={
                    iron_oct_delta: (0, 1, 2, 3),
                    bb3: (4, 5, 6, 7, 8, 9),
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

    This topology places a ditopic spacer between the vertices of the
    tetrahedron.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 3
        | 2-functional groups: 4 to 9

    Examples
    --------
    *Building Metal-Organic Tetrahedron*

    Many metal-organic cages are built using a process called
    subcomponent self-assembly, which is a complex chemical process
    that occurs in solution. In :class:`.Cage`, an example is
    provided of an alchemical approach using
    :class:`.M4L6TetrahedronSpacer` to
    build these types of cages. It is alchemical because the bonds
    formed during construction are not the same as the experimental
    reaction. Instead of forming bonds at the metal centre, we create
    bonds between disconnected ligands.

    The :class:`.M4L6TetrahedronSpacer` topology is provided for cases
    where the linkers cannot be disconnected symmetrically. However,
    in the case that the linker can be disconnected
    in a symmetrical fashion, the
    :class:`.M4L6Tetrahedron` topology can be used.

    """

    _non_linears = (
        NonLinearVertex(0, [0, 0, np.sqrt(6)/2]),
        NonLinearVertex(1, [-1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        NonLinearVertex(2, [1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        NonLinearVertex(3, [0, 2*np.sqrt(3)/3, -np.sqrt(6)/6]),
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
            vertices=(_non_linears[0], _non_linears[3]),
        ),
        LinearVertex.init_at_center(
            id=7,
            vertices=(_non_linears[1], _non_linears[2]),
        ),
        LinearVertex.init_at_center(
            id=8,
            vertices=(_non_linears[1], _non_linears[3]),
        ),
        LinearVertex.init_at_center(
            id=9,
            vertices=(_non_linears[2], _non_linears[3]),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[7]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[8]),
        Edge(6, _vertex_prototypes[2], _vertex_prototypes[5]),
        Edge(7, _vertex_prototypes[2], _vertex_prototypes[7]),
        Edge(8, _vertex_prototypes[2], _vertex_prototypes[9]),
        Edge(9, _vertex_prototypes[3], _vertex_prototypes[6]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[8]),
        Edge(11, _vertex_prototypes[3], _vertex_prototypes[9]),
    )

    _num_windows = 4
    _num_window_types = 1
