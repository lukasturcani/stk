"""
Cage
====

Organic
-------

.. toctree::
    :maxdepth: 2

    One Plus One <\
stk.molecular.topology_graphs.cage.three_plus_three.one_plus_one\
>
    Two Plus Two <\
stk.molecular.topology_graphs.cage.three_plus_three.two_plus_two\
>
    Two Plus Three <\
stk.molecular.topology_graphs.cage.two_plus_three.two_plus_three\
>
    Two Plus Four <\
stk.molecular.topology_graphs.cage.two_plus_four.two_plus_four\
>
    Three Plus Six <\
stk.molecular.topology_graphs.cage.two_plus_four.three_plus_six\
>
    Four Plus Four <\
stk.molecular.topology_graphs.cage.three_plus_three.four_plus_four\
>
    Four Plus Six <\
stk.molecular.topology_graphs.cage.two_plus_three.four_plus_six\
>
    Four Plus Six 2 <\
stk.molecular.topology_graphs.cage.two_plus_three.four_plus_six_2\
>
    Four Plus Eight <\
stk.molecular.topology_graphs.cage.two_plus_four.four_plus_eight\
>
    Five Plus Ten <\
stk.molecular.topology_graphs.cage.two_plus_four.five_plus_ten\
>
    Six Plus Eight <\
stk.molecular.topology_graphs.cage.three_plus_four.six_plus_eight\
>
    Six Plus Nine <\
stk.molecular.topology_graphs.cage.two_plus_three.six_plus_nine\
>
    Six Plus Twelve <\
stk.molecular.topology_graphs.cage.two_plus_four.six_plus_twelve\
>
    Eight Plus Twelve <\
stk.molecular.topology_graphs.cage.two_plus_three.eight_plus_twelve\
>
    Eight Plus Sixteen <\
stk.molecular.topology_graphs.cage.two_plus_four.eight_plus_sixteen\
>
    Ten Plus Twenty <\
stk.molecular.topology_graphs.cage.two_plus_four.ten_plus_twenty\
>
    Twelve Plus Thirty <\
stk.molecular.topology_graphs.cage.two_plus_five.twelve_plus_thirty\
>
    Twenty Plus Thirty <\
stk.molecular.topology_graphs.cage.two_plus_three.twenty_plus_thirty\
>

Metal-Organic
-------------

.. toctree::
    :maxdepth: 2

    M2L4 Lantern <\
stk.molecular.topology_graphs.cage.metal_topologies.m2l4_lantern\
>
    M3L3 Triangle <\
stk.molecular.topology_graphs.cage.metal_topologies.m3l3_triangle\
>
    M3L6 <\
stk.molecular.topology_graphs.cage.metal_topologies.m3l6\
>
    M4L4 Square <\
stk.molecular.topology_graphs.cage.metal_topologies.m4l4_square\
>
    M4L4 Tetrahedron <\
stk.molecular.topology_graphs.cage.metal_topologies.m4l4_tetrahedron\
>
    M4L6 Tetrahedron Spacer <\
stk.molecular.topology_graphs.cage.metal_topologies\
.m4l6_tetrahedron_spacer\
>
    M4L6 Tetrahedron <\
stk.molecular.topology_graphs.cage.metal_topologies.m4l6_tetrahedron\
>
    M4L8 <\
stk.molecular.topology_graphs.cage.metal_topologies.m4l8\
>
    M6L2L3 Prism <\
stk.molecular.topology_graphs.cage.metal_topologies.m6l2l3_prism\
>
    M6L12 Cube <\
stk.molecular.topology_graphs.cage.metal_topologies.m6l12_cube\
>
    M8L6 Cube <\
stk.molecular.topology_graphs.cage.metal_topologies.m8l6_cube\
>
    M12L24 <\
stk.molecular.topology_graphs.cage.metal_topologies.m12l24\
>
    M24L48 <\
stk.molecular.topology_graphs.cage.metal_topologies.m24l48\
>

"""

from collections import Counter, defaultdict
from functools import partial

from ...reactions import GenericReactionFactory
from ..topology_graph import NullOptimizer, TopologyGraph
from .cage_construction_state import _CageConstructionState
from .vertices import UnaligningVertex


class UnoccupiedVertexError(Exception):
    """
    When a cage vertex is not occupied by a building block.

    """

    pass


class OverlyOccupiedVertexError(Exception):
    """
    When a cage vertex is occupied by more than one building block.

    """

    pass


class Cage(TopologyGraph):
    """
    Represents a cage topology graph.

    Notes
    -----
    Cage topologies are added by creating a subclass, which defines the
    :attr:`_vertex_prototypes` and :attr:`_edge_prototypes` class
    attributes.

    .. _cage-topology-graph-examples:

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in :mod:`~.cage.cage`,
    can serve as good examples.

    *Basic Construction*

    :class:`.Cage` instances can be made by providing the building
    block molecules only (using :class:`.FourPlusSix` as an example)

    .. testcode:: basic-construction

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
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
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )

    *Suggested Optimization*

    For :class:`.Cage` topologies, it is recommend to use the
    :class:`.MCHammer` optimizer.
    However, for cages formed from highly unsymmetrical building
    blocks, it is recommend to use the simplified
    :class:`.Collapser` optimizer.


    .. testcode:: suggested-optimization

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(bb1, bb2),
                optimizer=stk.MCHammer(),
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(bb1, bb2),
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
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )

    *Structural Isomer Construction*

    Different structural isomers of cages can be made by using the
    `vertex_alignments` optional parameter

    .. testcode:: structural-isomer-construction

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(bb1, bb2),
                vertex_alignments={0: 1, 1: 1, 2: 2},
            ),
        )

    The parameter maps the id of a vertex to a number
    between 0 (inclusive) and the number of edges the vertex is
    connected to (exclusive). So a vertex connected to three edges
    can be mapped to ``0``, ``1`` or ``2``.

    By changing which edge each vertex is aligned with, a different
    structural isomer of the cage can be formed.

    *Multi-Building Block Cage Construction*

    You can also build cages with multiple building blocks, but,
    if you have multiple building blocks with the same number
    of functional groups, you have to assign each building block to the
    vertex you want to place it on

    .. testcode:: multi-building-block-cage-construction

        import stk

        bb1 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='O=CC(Cl)(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb3 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb4 = stk.BuildingBlock(
            smiles='NCC(Cl)N',
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        bb5 = stk.BuildingBlock('NCCCCN', [stk.PrimaryAminoFactory()])

        cage1 = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks={
                    bb1: range(2),
                    bb2: (2, 3),
                    bb3: 4,
                    bb4: 5,
                    bb5: range(6, 10),
                },
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='O=CC(Cl)(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb3 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb4 = stk.BuildingBlock(
            smiles='NCC(Cl)N',
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        bb5 = stk.BuildingBlock('NCCCCN', [stk.PrimaryAminoFactory()])

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks={
                    bb1: range(2),
                    bb2: (2, 3),
                    bb3: 4,
                    bb4: 5,
                    bb5: range(6, 10),
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
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )


    You can combine this with the `vertex_alignments` parameter

    .. testcode:: multi-building-block-cage-construction

        cage2 = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks={
                    bb1: range(2),
                    bb2: (2, 3),
                    bb3: 4,
                    bb4: 5,
                    bb5: range(6, 10),
                },
                vertex_alignments={5: 1},
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='O=CC(Cl)(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb3 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb4 = stk.BuildingBlock(
            smiles='NCC(Cl)N',
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        bb5 = stk.BuildingBlock('NCCCCN', [stk.PrimaryAminoFactory()])

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks={
                    bb1: range(2),
                    bb2: (2, 3),
                    bb3: 4,
                    bb4: 5,
                    bb5: range(6, 10),
                },
                vertex_alignments={5: 1},
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
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )

    *Metal-Organic Cage Construction*

    A series of common metal-organic cage topologies are provided and
    can be constructed in the same way as other :class:`.Cage`
    instances using metal atoms and :class:`DativeReactionFactory`
    instances to produce metal-ligand bonds. Each metal topology has
    specific vertices reserved for the metal atoms or complexes,
    which are listed in their documentation.

    .. testcode:: metal-organic-cage-construction

        import stk

        # Produce a Pd+2 atom with 4 functional groups.
        palladium_atom = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=[[0., 0., 0.]],
        )

        # Build a building block with two functional groups using
        # the SmartsFunctionalGroupFactory.
        bb1 = stk.BuildingBlock(
            smiles=(
                'C1=NC=CC(C2=CC=CC(C3=C'
                'C=NC=C3)=C2)=C1'
            ),
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        cage1 = stk.ConstructedMolecule(
            stk.cage.M2L4Lantern(
                building_blocks=(palladium_atom, bb1),
                # Ensure that bonds between the GenericFunctionalGroups
                # of the ligand and the SingleAtom functional groups
                # of the metal are dative.
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            frozenset({
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom,
                            }): 9,
                        },
                    ),
                ),
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        palladium_atom = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=[[0., 0., 0.]],
        )

        bb1 = stk.BuildingBlock(
            smiles=(
                'C1=NC=CC(C2=CC=CC(C3=C'
                'C=NC=C3)=C2)=C1'
            ),
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        cage = stk.ConstructedMolecule(
            stk.cage.M2L4Lantern(
                building_blocks=(palladium_atom, bb1),
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            # Use bond order of 1 here so that the
                            # rendering does not show a bond order
                            # of 9.
                            frozenset({
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom,
                            }): 1,
                        },
                    ),
                ),
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
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )

    *Controlling Metal-Complex Stereochemistry*

    When building metal-organic cages from octahedral metals, i.e.
    Fe(II), the stereochemistry of the metal centre can be important.
    Maintaining that stereochemistry around specific metal centres
    during :class:`.Cage` construction is difficult, so an
    alternative route to these types of structures can be taken.
    Firstly, you would construct a :class:`.MetalComplex` instance
    with the appropriate stereochemistry and dummy reactive groups
    (bromine in the following example)

    .. testcode:: controlling-metal-complex-stereochemistry

        import stk

        # Produce a Fe+2 atom with 6 functional groups.
        iron_atom = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        # Define coordinating ligand with dummy bromine groups and
        # metal coordinating functional groups.
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

        # Build iron complex with delta stereochemistry.
        iron_oct_delta = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralDelta(
                metals=iron_atom,
                ligands=bb2,
            ),
        )

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
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    complex.get_atoms(),
                    complex.get_position_matrix(),
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
                ) for bond in complex.get_bonds()
            ),
        )

    Then the metal complexes can be placed on the appropriate
    :class:`.Cage` topology to produce a structure with the desired
    stereochemistry at all metal centres.

    .. testcode:: controlling-metal-complex-stereochemistry

        # Assign Bromo functional groups to the metal complex.
        iron_oct_delta = stk.BuildingBlock.init_from_molecule(
            molecule=iron_oct_delta,
            functional_groups=[stk.BromoFactory()],
        )

        # Define spacer building block.
        bb3 = stk.BuildingBlock(
            smiles=(
                'C1=CC(C2=CC=C(Br)C=C2)=C'
                'C=C1Br'
            ),
            functional_groups=[stk.BromoFactory()],
        )

        # Build an M4L6 Tetrahedron with a spacer.
        cage2 = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L6TetrahedronSpacer(
                building_blocks=(
                    iron_oct_delta,
                    bb3,
                ),
            ),
        )

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

        iron_oct_delta = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralDelta(
                metals=iron_atom,
                ligands=bb2,
            ),
        )

        iron_oct_delta = stk.BuildingBlock.init_from_molecule(
            molecule=iron_oct_delta,
            functional_groups=[stk.BromoFactory()],
        )

        bb3 = stk.BuildingBlock(
            smiles=(
                'C1=CC(C2=CC=C(Br)C=C2)=C'
                'C=C1Br'
            ),
            functional_groups=[stk.BromoFactory()],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L6TetrahedronSpacer(
                building_blocks=(
                    iron_oct_delta,
                    bb3,
                ),
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

    *Aligning Metal Complex Building Blocks*

    When building metal-organic cages from metal complex
    building blocks, it is common that
    the metal complex :class:`.BuildingBlock` will have
    multiple functional groups, but that those functional groups
    are overlapping. This means that some of its atoms appear in
    multiple functional groups. A difficulty arises when the
    atom shared between the functional groups is a *placer* atom.

    *Placer* atoms are used to align building blocks, so that
    they have an appropriate orientation in the final topology.
    If there is only one *placer* atom, no alignment can be made,
    as no vector running between *placer* atoms can be defined,
    and used for the alignment of the :class:`.BuildingBlock`.

    By default, :mod:`stk` may create overlapping functional
    groups, which may lead to a lack of an appropriate number
    of *placer* atoms, leading to a :class:`.BuildingBlock`
    being unaligned. However, the user can manually set the
    *placer* atoms of functional groups, so that not all of the
    *placer* atoms appear in multiple functional groups, which
    leads to proper alignment.

    First we build a metal complex

    .. testcode:: aligning-metal-complex-building-blocks

        import stk

        metal_atom = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=[[0., 0., 0.]],
        )

        ligand = stk.BuildingBlock(
            smiles='NCCN',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#7]~[#6]',
                    bonders=(0, ),
                    deleters=(),
                ),
            ],
        )

        metal_complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
                metals=metal_atom,
                ligands=ligand,
            ),
        )

    Next, we convert the metal complex into a :class:`.BuildingBlock`,
    taking care to define functional groups which do not have
    overlapping *placer* atoms

    .. testcode:: aligning-metal-complex-building-blocks

        metal_complex = stk.BuildingBlock.init_from_molecule(
            molecule=metal_complex,
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[Pd]~[#7]',
                    bonders=(0, ),
                    deleters=(),
                    # The nitrogen atom will be different
                    # for each functional group.
                    placers=(0, 1),
                ),
            ],
        )

    We load in the organic linker of the cage as normal

    .. testcode:: aligning-metal-complex-building-blocks

        linker = stk.BuildingBlock(
            smiles='C1=NC=CC(C2=CC=NC=C2)=C1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

    And finally, we build the cage with a
    :class:`DativeReactionFactory` instance to produce dative bonds.

    .. testcode:: aligning-metal-complex-building-blocks

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L4Square(
                corners=metal_complex,
                linkers=linker,
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            frozenset({
                                stk.GenericFunctionalGroup,
                                stk.GenericFunctionalGroup,
                            }): 9,
                        },
                    ),
                ),
            ),
        )

    """

    def __init_subclass__(cls, **kwargs):
        cls._vertex_degrees = Counter(
            vertex_id
            for edge in cls._edge_prototypes
            for vertex_id in edge.get_vertex_ids()
        )
        cls._vertices_of_degree = defaultdict(set)
        for vertex_id, degree in cls._vertex_degrees.items():
            cls._vertices_of_degree[degree].add(vertex_id)

    def __init__(
        self,
        building_blocks,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
        optimizer=NullOptimizer(),
    ):
        """
        Initialize a :class:`.Cage`.

        Parameters
        ----------
        building_blocks : :class:`iterable` or :class:`dict`
            Can be a :class:`iterable` of :class:`.BuildingBlock`
            instances, which should be placed on the topology graph.

            Can also be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on. A :class:`dict` is
            required when there are multiple building blocks with the
            same number of functional groups, because in this case
            the desired placement is ambiguous.

        vertex_alignments : :class:`dict`, optional
            A mapping from the id of a :class:`.Vertex`
            to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is referred to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        optimizer : :class:`.Optimizer`, optional
            Used to optimize the structure of the constructed
            molecule.

        Raises
        ------
        :class:`AssertionError`
            If the any building block does not have a
            valid number of functional groups.

        :class:`ValueError`
            If the there are multiple building blocks with the
            same number of functional_groups in `building_blocks`,
            and they are not explicitly assigned to vertices. The
            desired placement of building blocks is ambiguous in
            this case.

        :class:`~.cage.UnoccupiedVertexError`
            If a vertex of the cage topology graph does not have a
            building block placed on it.

        :class:`~.cage.OverlyOccupiedVertexError`
            If a vertex of the cage topology graph has more than one
            building block placed on it.

        """

        building_block_vertices = self._normalize_building_blocks(
            building_blocks=building_blocks,
        )
        self._vertex_alignments = (
            dict(vertex_alignments)
            if vertex_alignments is not None
            else {}
        )
        building_block_vertices = self._with_unaligning_vertices(
            building_block_vertices=building_block_vertices,
        )
        building_block_vertices = self._assign_aligners(
            building_block_vertices=building_block_vertices,
            vertex_alignments=self._vertex_alignments,
        )
        self._check_building_block_vertices(building_block_vertices)
        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=self._edge_prototypes,
            reaction_factory=reaction_factory,
            construction_stages=tuple(
                partial(self._has_degree, degree)
                for degree
                in sorted(self._vertices_of_degree, reverse=True)
            ),
            num_processes=num_processes,
            optimizer=optimizer,
            edge_groups=None,
        )

    @classmethod
    def _normalize_building_blocks(cls, building_blocks):
        # Use tuple here because it prints nicely.
        allowed_degrees = tuple(cls._vertices_of_degree.keys())
        if isinstance(building_blocks, dict):
            for building_block in building_blocks:
                assert (
                    building_block.get_num_functional_groups()
                    in cls._vertices_of_degree.keys()
                ), (
                    'The number of functional groups in '
                    f'{building_block} needs to be one of '
                    f'{allowed_degrees}, but is '
                    'currently '
                    f'{building_block.get_num_functional_groups()}.'
                )
            return {
                building_block: cls._get_vertices(ids)
                for building_block, ids in building_blocks.items()
            }

        else:
            return cls._get_building_block_vertices(
                building_blocks=building_blocks,
            )

    @staticmethod
    def _with_unaligning_vertices(building_block_vertices):
        clone = dict(building_block_vertices)
        for building_block, vertices in clone.items():
            # Building blocks with 1 placer, cannot be aligned and
            # must therefore use an UnaligningVertex.
            if building_block.get_num_placers() == 1:
                clone[building_block] = tuple(
                    UnaligningVertex(
                        id=v.get_id(),
                        position=v.get_position(),
                        use_neighbor_placement=(
                            v.use_neighbor_placement()
                        ),
                        aligner_edge=v.get_aligner_edge(),
                    )
                    for v in vertices
                )

        return clone

    @classmethod
    def _assign_aligners(
        cls,
        building_block_vertices,
        vertex_alignments,
    ):
        def with_aligner(vertex):
            return vertex.with_aligner_edge(
                aligner_edge=vertex_alignments.get(vertex.get_id(), 0),
            )

        return {
            building_block: tuple(map(with_aligner, vertices))
            for building_block, vertices
            in building_block_vertices.items()
        }

    @classmethod
    def _check_building_block_vertices(cls, building_block_vertices):
        unassigned_ids = set(
            vertex.get_id() for vertex in cls._vertex_prototypes
        )
        assigned_ids = set()
        vertices = (
            vertex
            for vertices_ in building_block_vertices.values()
            for vertex in vertices_
        )
        for vertex in vertices:
            if vertex.get_id() in assigned_ids:
                raise OverlyOccupiedVertexError(
                    f'Vertex {vertex.get_id()} has multiple building '
                    'blocks placed on it.'
                )
            assigned_ids.add(vertex.get_id())
            unassigned_ids.remove(vertex.get_id())

        if unassigned_ids:
            raise UnoccupiedVertexError(
                'The following vertices are unoccupied '
                f'{unassigned_ids}.'
            )

    def clone(self):
        clone = super().clone()
        clone._vertex_alignments = dict(self._vertex_alignments)
        return clone

    @classmethod
    def _get_vertices(cls, vertex_ids):
        """
        Yield vertex prototypes.

        Parameters
        ----------
        vertex_ids : :class:`iterable` of :class:`int`
            The ids of the vertices to yield.

        Yields
        ------
        :class:`.Vertex`
            A vertex prototype of the topology graph.

        """

        if isinstance(vertex_ids, int):
            vertex_ids = (vertex_ids, )

        for vertex_id in vertex_ids:
            yield cls._vertex_prototypes[vertex_id]

    def _has_degree(self, degree, vertex):
        """
        Check if `vertex` has a degree of `degree`.

        Parameters
        ----------
        degree : :class:`int`
            The degree in question.

        vertex : :class:`.Vertex`
            The vertex in question.

        Returns
        -------
        :class:`bool`
            ``True`` if `vertex` has a degree of `degree`.

        """

        return vertex.get_id() in self._vertices_of_degree[degree]

    @classmethod
    def _get_building_block_vertices(cls, building_blocks):
        """
        Map building blocks to the vertices of the graph.

        Parameters
        ----------
        building_blocks : :class:`iterable` of :class:`.BuildingBlock`
            The building blocks which need to be mapped to vertices.

        Returns
        -------
        :class:`dict`
            Maps each building block in `building_blocks` to a
            :class:`list` of :class:`.Vertex` instances it should be
            placed on.

        Raises
        ------
        :class:`AssertionError`
            If the any building block does not have a
            valid number of functional groups.

        :class:`ValueError`
            If there are multiple building blocks with the same number
            of functional groups.

        """

        # Use tuple here because it prints nicely.
        allowed_degrees = tuple(cls._vertices_of_degree.keys())

        building_blocks_by_degree = {}
        for building_block in building_blocks:
            num_fgs = building_block.get_num_functional_groups()
            assert (
                num_fgs in cls._vertices_of_degree.keys()
            ), (
                'The number of functional groups in '
                f'{building_block} needs to be one of '
                f'{allowed_degrees}, but is '
                'currently '
                f'{building_block.get_num_functional_groups()}.'
            )

            if num_fgs in building_blocks_by_degree:
                raise ValueError(
                    'If there are multiple building blocks with the '
                    'same number of functional groups, '
                    'building_block_vertices must be set explicitly.'
                )
            building_blocks_by_degree[num_fgs] = building_block

        building_block_vertices = {}
        for vertex in cls._vertex_prototypes:
            vertex_degree = cls._vertex_degrees[vertex.get_id()]
            building_block = building_blocks_by_degree[vertex_degree]
            building_block_vertices[building_block] = (
                building_block_vertices.get(building_block, [])
            )
            building_block_vertices[building_block].append(vertex)
        return building_block_vertices

    def _get_scale(self, building_block_vertices):
        return max(
            bb.get_maximum_diameter()
            for bb in building_block_vertices
        )

    def _get_construction_state(self):
        return _CageConstructionState(
            building_block_vertices=self._building_block_vertices,
            edges=self._edges,
            num_placement_stages=self._implementation.get_num_stages(),
            vertex_degrees=self._vertex_degrees,
        )

    def __repr__(self):
        vertex_alignments = (
            f'vertex_alignments={self._vertex_alignments}'
            if self._vertex_alignments
            else ''
        )
        return f'cage.{self.__class__.__name__}({vertex_alignments})'
