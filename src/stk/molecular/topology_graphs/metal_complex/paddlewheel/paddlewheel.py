"""
Paddlewheel
===========

"""

from ...topology_graph import Edge
from ..metal_complex import MetalComplex
from ..vertices import BiDentateLigandVertex, MetalVertex


class Paddlewheel(MetalComplex):
    """
    Represents a metal complex topology graph.

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='[Cu+2]',
            functional_groups=(
                stk.SingleAtom(stk.Cu(0, charge=2))
                for i in range(4)
            ),
            position_matrix=([0, 0, 0], ),
        )
        bb2 = stk.BuildingBlock(
            smiles='O=C(O)c1ccc(Br)cc1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#8]~[#1]',
                    bonders=(1, ),
                    deleters=(2, ),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#8X1]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )
        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.Paddlewheel(
                metals=bb1,
                ligands=bb2,
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

    Metal building blocks with at least four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology graph.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: (0, 1)
        | ligands: (0, 1, 2, 3)

    See :class:`.MetalComplex` for more details and examples.

    """

    _metal_vertex_prototypes = (
        MetalVertex(0, [0, 1, 0]),
        MetalVertex(1, [0, -1, 0]),
    )
    _ligand_vertex_prototypes = (
        BiDentateLigandVertex(2, [2, 0, 0]),
        BiDentateLigandVertex(3, [0, 0, 2]),
        BiDentateLigandVertex(4, [-2, 0, 0]),
        BiDentateLigandVertex(5, [0, 0, -2]),
    )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0.1, 0.5, 0],
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[1],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0.1, -0.5, 0],
        ),

        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
            position=[0, 0.5, 0.1],
        ),
        Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[1],
            vertex2=_ligand_vertex_prototypes[1],
            position=[0, -0.5, 0.1],
        ),

        Edge(
            id=4,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=[-0.1, 0.5, 0],
        ),
        Edge(
            id=5,
            vertex1=_metal_vertex_prototypes[1],
            vertex2=_ligand_vertex_prototypes[2],
            position=[-0.1, -0.5, 0],
        ),

        Edge(
            id=6,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[3],
            position=[0, 0.5, -0.1],
        ),
        Edge(
            id=7,
            vertex1=_metal_vertex_prototypes[1],
            vertex2=_ligand_vertex_prototypes[3],
            position=[0, -0.5, -0.1],
        ),
    )
