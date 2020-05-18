"""
M4L4 Square
===========

"""

from ..cage import Cage
from ..vertices import _LinearCageVertex
from ...topology_graph import Edge
from ...reactions import GenericReactionFactory


class M4L4Square(Cage):
    """
    Represents a cage topology graph.

    Both `corner` and `linker` vertices require building blocks with
    two functional groups for this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | corners: (0, 1, 2, 3)
        | linkers: (4, 5, 6, 7)

    Examples
    --------
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
    being unalinged. However, the user can manually set the 
    *placer* atoms of functional groups, so that the *placer*
    atoms do not appear in multiple functional groups, which
    leads to proper alignment.
    
    First we build a metal complex

    .. code-block:: python

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
            ]
        )

        metal_complex = stk.ConstructedMolecule(
            stk.metal_complex.CisProtectedSquarePlanar(
                metals=metal_atom,
                ligands=ligand,
            )
        )
        
    Next, we convert the metal complex into a :class:`.BuildingBlock`,
    taking care to define functional groups which do not have
    overlapping *placer* atoms

    .. code-block:: python

        metal_complex = stk.BuildingBlock.init_from_molecule(
            molecule=metal_complex,
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[Pd]~[#7]',
                    bonders=(0, ),
                    deleters=(),
                    # The nitrogen atom will be different
                    # for each functional group.
                    placers=(1, ),
                ),
            ]
        )

    We load in the organic linker of the cage as normal

    .. code-block:: python

        linker = stk.BuildingBlock(
            smiles='C1=NC=CC(C2=CC=NC=C2)=C1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ]
        )

    And finally, we build the cage with a
    :class:`DativeReactionFactory` instance to produce dative bonds.

    .. code-block:: python

        cage = stk.ConstructedMolecule(
            stk.cage.M4L4Square(
                corners=metal_complex,
                linkers=linker,
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            frozenset({
                                stk.GenericFunctionalGroup,
                                stk.GenericFunctionalGroup,
                            }): 9
                        }
                    )
                ),
            )
        )

    See :class:`.Cage` for more details and examples.

    """

    def __init__(
        self,
        corners,
        linkers,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):
        """
        Initialize a :class:`.M4L4Square`.

        Parameters
        ----------
        corners : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed on all corner vertices on the topology
            graph.

        linkers : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed on all linker vertices on the topology
            graph.

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

        """

        if isinstance(corners, dict):
            building_blocks = {
                corner: tuple(self._get_vertices(ids))
                for corner, ids in corners.items()
            }
        else:
            ids = (0, 1, 2, 3)
            building_blocks = {
                corners: tuple(self._get_vertices(ids))
            }

        if isinstance(linkers, dict):
            linkers_dict = {
                linker: tuple(self._get_vertices(ids))
                for linker, ids in corners.items()
            }
        else:
            ids = (4, 5, 6, 7)
            linkers_dict = {
                linkers: tuple(self._get_vertices(ids))
            }

        building_blocks.update(
            (building_block, vertices)
            for building_block, vertices in linkers_dict.items()
        )

        super().__init__(
            building_blocks,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
        )

    _vertex_prototypes = (
        _LinearCageVertex(0, [1, 1, 0]),
        _LinearCageVertex(1, [1, -1, 0]),
        _LinearCageVertex(2, [-1, -1, 0]),
        _LinearCageVertex(3, [-1, 1, 0]),

        _LinearCageVertex(4, [1, 0, 0], False),
        _LinearCageVertex(5, [0, -1, 0], False),
        _LinearCageVertex(6, [-1, 0, 0], False),
        _LinearCageVertex(7, [0, 1, 0], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[4]),

        Edge(2, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(3, _vertex_prototypes[2], _vertex_prototypes[5]),

        Edge(4, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(5, _vertex_prototypes[3], _vertex_prototypes[6]),

        Edge(6, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(7, _vertex_prototypes[0], _vertex_prototypes[7]),
    )

    _num_windows = 1
    _num_window_types = 1
