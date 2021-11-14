"""
M4L6 Tetrahedron
================

"""

import numpy as np

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import NonLinearVertex


class M4L6Tetrahedron(Cage):
    """
    Represents a cage topology graph.

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
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

        iron_atom = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        iron_oct_delta = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralDelta(
                metals=iron_atom,
                ligands=bb1,
            ),
        )

        iron_oct_delta = stk.BuildingBlock.init_from_molecule(
            molecule=iron_oct_delta,
            functional_groups=[stk.BromoFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L6Tetrahedron(
                building_blocks=(iron_oct_delta, ),
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

        bb1 = stk.BuildingBlock(
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

        iron_atom = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        iron_oct_delta = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralDelta(
                metals=iron_atom,
                ligands=bb1,
                optimizer=stk.MCHammer(),
            ),
        )

        iron_oct_delta = stk.BuildingBlock.init_from_molecule(
            molecule=iron_oct_delta,
            functional_groups=[stk.BromoFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L6Tetrahedron(
                building_blocks=(iron_oct_delta, ),
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

    This topology directly connects the vertices of the tetrahedron.

    Building blocks with three functional groups are required for this
    topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 3

    Examples
    --------
    *Building a Metal-Organic Tetrahedron*

    Many metal-organic cages are built using a process called
    subcomponent self-assembly, which is a complex chemical process
    that occurs in solution. Here, we provide an example of an
    alchemical approach to building these types of cages. It is
    alchemical because the bonds formed during construction are not
    the same as the experimental reaction. Instead of forming bonds at
    the metal centre, we create bonds between disconnected ligands.
    Firstly, we define the disconnected linker molecule, where we have
    added a bromine atom at the disconnection to perform the cage
    reaction.

    .. testcode:: building-a-metal-organic-tetrahedron

        import stk

        # Define coordinating ligand with dummy bromine groups and
        # metal coordinating functional groups.
        bb1 = stk.BuildingBlock(
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

    We then build the desired metal complex with this linker. This
    process completes all metal-ligand reactions performed during the
    subcomponent self-assembly process.

    .. testcode:: building-a-metal-organic-tetrahedron

        # Produce a Fe+2 atom with 6 functional groups.
        iron_atom = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        # Build iron complex with delta stereochemistry.
        iron_oct_delta = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralDelta(
                metals=iron_atom,
                ligands=bb1,
            ),
        )

    We then redefine the metal complex building block based on its
    bromine functional groups, which becomes the metal-based
    building block of any cage formed by this process.

    .. testcode:: building-a-metal-organic-tetrahedron

        # Assign Bromo functional groups to the metal complex.
        iron_oct_delta = stk.BuildingBlock.init_from_molecule(
            molecule=iron_oct_delta,
            functional_groups=[stk.BromoFactory()],
        )

    Finally, we build the :class:`M4L6Tetrahedron` cage using this
    building block.

    .. testcode:: building-a-metal-organic-tetrahedron

        # Build an M4L6 Tetrahedron.
        cage2 = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L6Tetrahedron(
                building_blocks=(iron_oct_delta, ),
            ),
        )

    Importantly, in the case that the linker cannot be disconnected
    in a symmetrical fashion, we have provided the
    :class:`M4L6TetrahedronSpacer` topology, which has a spacer
    vertex between the metal vertices on the tetrahedron. See
    :class:`.M4L6TetrahedronSpacer` for an example of its usage.

    See :class:`.Cage` for more details and examples.

    """

    _x = 1/(2*np.sqrt(2))
    _y = 0.5
    _vertex_prototypes = (
        NonLinearVertex(0, [_y, 0, -_x]),
        NonLinearVertex(1, [-_y, 0, -_x]),
        NonLinearVertex(2, [0, _y, _x]),
        NonLinearVertex(3, [0, -_y, _x]),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[1]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(5, _vertex_prototypes[2], _vertex_prototypes[3]),
    )

    _num_windows = 4
    _num_window_types = 1
