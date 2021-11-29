"""
Octahedral
==========

"""

from ...topology_graph import Edge
from ..metal_complex import MetalComplex
from ..vertices import MetalVertex, MonoDentateLigandVertex


class Octahedral(MetalComplex):
    """
    Represents a metal complex topology graph.

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )
        bb2 = stk.BuildingBlock(
            smiles='C1=CC=NC=C1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )
        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.Octahedral(
                metals=bb1,
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

    Metal building blocks with at least six functional groups are
    required for this topology.

    Ligand building blocks with one functional group are required for
    this topology graph.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: (0, )
        | ligands: (0, 1, 2, 3, 4, 5)

    See :class:`.MetalComplex` for more details and examples.

    """

    _metal_vertex_prototypes = (
        MetalVertex(0, [0, 0, 0]),
    )
    _ligand_vertex_prototypes = (
        MonoDentateLigandVertex(1, [2.5, 0, 0]),
        MonoDentateLigandVertex(2, [0, 2.5, 0]),
        MonoDentateLigandVertex(3, [0, 0, 2.5]),
        MonoDentateLigandVertex(4, [-2.5, 0, 0]),
        MonoDentateLigandVertex(5, [0, -2.5, 0]),
        MonoDentateLigandVertex(6, [0, 0, -2.5]),
    )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
        ),
        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
        ),
        Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[3],
        ),
        Edge(
            id=4,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[4],
        ),
        Edge(
            id=5,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[5],
        ),
    )
