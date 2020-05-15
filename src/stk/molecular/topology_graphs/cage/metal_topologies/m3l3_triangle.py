"""
M3L3 Triangle
=============

"""

import numpy as np

from ..cage import Cage
from ..vertices import _LinearCageVertex
from ...topology_graph import Edge
from ...reactions import GenericReactionFactory


class M3L3Triangle(Cage):
    """
    Represents a cage topology graph.

    Both `corner` and `linker` vertices require building blocks with
    two functional groups for this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | corners: (0, 1, 2)
        | linkers: (3, 4, 5)

    Examples
    --------

    *Advanced Construction*

    When building metal-organic structures, it is common that
    functional groups, defined by metal atoms, will be overlapping.
    Therefore, to build :class:`.M3L3Triangle` instances with the
    desired building block alignment, corner and linker building blocks
    must have functional groups with additional `placer` atoms. First
    we build the metal atom (Pd(II) with four functional groups).


    .. code-block:: python

        import stk

        metal_atom = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=([0, 0, 0], ),
        )

    Next we build the corner unit, which is a cis-protected square
    planar metal complex.

    .. code-block:: python

        bb1 = stk.BuildingBlock(
            smiles='NCCN',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#7]~[#6]',
                    bonders=(0, ),
                    deleters=(),
                ),
            ]
        )

        corner_bb = stk.ConstructedMolecule(
            stk.metal_complex.CisProtectedSquarePlanar(
                metals=metal_atom,
                ligands=bb1,
            )
        )

    We then define a new :class:`.BuildingBlock` with functional groups
    that have two placer atoms to enforce alignment of the corner unit.

    .. code-block:: python

        corner_bb = stk.BuildingBlock.init_from_molecule(
            molecule=corner_bb,
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[Pd]~[#7]',
                    bonders=(0, ),
                    deleters=(),
                    placers=(0, 1),
                ),
            ]
        )

    We also load in the organic linker as normal.

    .. code-block:: python

        linker_bb = stk.BuildingBlock(
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
    The intializer for this topology is different to other
    :class:`.Cage` instances because both building block types have
    the same number of functional groups.

    .. code-block:: python

        cage = stk.ConstructedMolecule(
            stk.cage.M3L3Triangle(
                corners=corner_bb,
                linkers=linker_bb,
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            frozenset({
                                stk.GenericFunctionalGroup,
                                stk.GenericFunctionalGroup
                            }): 9
                        }
                    )
                )
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
        Initialize a :class:`.M3L3Triangle`.

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
            ids = (0, 1, 2)
            building_blocks = {
                corners: tuple(self._get_vertices(ids))
            }

        if isinstance(linkers, dict):
            linkers_dict = {
                linker: tuple(self._get_vertices(ids))
                for linker, ids in corners.items()
            }
        else:
            ids = (3, 4, 5)
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

    _x = 2*np.sqrt(3)/4
    _y = 2
    _vertex_prototypes = (
        _LinearCageVertex(0, [0, _x, 0]),
        _LinearCageVertex(1, [_y/2, -_x, 0]),
        _LinearCageVertex(2, [-_y/2, -_x, 0]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _LinearCageVertex.init_at_center(
            id=3,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=4,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=5,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[0]),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[3]),

        Edge(2, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[2], _vertex_prototypes[4]),

        Edge(4, _vertex_prototypes[2], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[0], _vertex_prototypes[5]),
    )

    _num_windows = 1
    _num_window_types = 1
