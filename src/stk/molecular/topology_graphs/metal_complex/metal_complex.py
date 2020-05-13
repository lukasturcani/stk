"""
Metal Complex
=============

.. toctree::
    :maxdepth: 2

    Paddlewheel <\
    stk.molecular.topology_graphs.metal_complex.paddlewheel.paddlewheel\
>
    Porphyrin <\
stk.molecular.topology_graphs.metal_complex.porphyrin.porphyrin\
>
    Octahedral <\
stk.molecular.topology_graphs.metal_complex.octahedral.octahedral\
>
    Octahedral Lambda <\
stk.molecular.topology_graphs.metal_complex.octahedral.octahedral_lambda\
>
    Octahedral Delta <\
stk.molecular.topology_graphs.metal_complex.octahedral.octahedral_delta\
>
    Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar.square_planar\
>
    Bidentate Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar.bidentate_square_planar\
>
    Cis Protected Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar.cis_protected_square_planar\
>

"""

from ..topology_graph import TopologyGraph
from ...reactions import GenericReactionFactory


class MetalComplex(TopologyGraph):
    """
    Represents a metal complex topology graph.

    Notes
    -----
    *Subclass Implementation*

    Each subclass needs to define the iterable attributes,
    :attr:`_metal_vertex_prototypes` and
    :attr:`_ligand_vertex_prototypes`, which are iterables of
    :class:`.Vertex` instances and define the positions for
    metal-based and ligand-based vertices, respectively.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in
    :mod:`~.metal_complex.metal_complex`, can serve as good examples.

    *Basic Construction*

    For most :class:`.MetalComplex` topology graphs, we first
    need to define a metal :class:`.BuildingBlock`, consisting of
    1 atom and multiple functional groups

    .. code-block:: python

        import stk
        import numpy as np

        metal = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=np.array([[0, 0, 0]]),
        )

    We also need to define an organic ligand :class:`.BuildingBlock`

    .. code-block:: python

        # Define an organic linker with two functional groups.
        bidentate = stk.BuildingBlock(
            smiles='C=NC/C=N/Br',
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
            ]
        )

    Finally, we can create the :class:`.MetalComplex`.

    .. code-block:: python

        complex = stk.ConstructedMolecule(
            stk.metal_complex.OctahedralLambda(
                metals=metal,
                ligands=bidentate,
            )
        )

    *Construction with Multiple Metals & Ligands*

    When multiple metals or ligands are used, the `metals` and
    `ligands` parameters accept values of type :class:`dict`, which
    specify the exact vertex each metal or ligand needs to be placed
    on. However, if each ligand is has a different number of
    functional groups, they can be provided together in a
    :class:`tuple`.

    Note that the valid vertex identifiers depend on the exact
    metal complex you are using. These are detailed in the docstring
    for that specific metal vertex topology graph.

    *Leaving Unsubstitued Sites*

    Sometimes you may want to build a metal comlex with open metal
    sites. For example, here we show how to build a square planar
    palladium(II) complex with two open metal sites.

    .. code-block:: python

        import stk

        pd = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=np.array([[0, 0, 0]]),
        )

        # Define a bidentate ligand with two functional groups.
        bidentate_ligand = stk.BuildingBlock(
            smiles='NCCN',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#7]~[#6]',
                    bonders=(0, ),
                    deleters=(),
                ),
            ]
        )

        # Construct a cis-protected square planar metal complex.
        complex = stk.ConstructedMolecule(
            stk.metal_complex.CisProtectedSquarePlanar(
                metals=pd_metal,
                ligands=bidentate_ligand,
            )
        )

    """

    def __init__(
        self,
        metals,
        ligands,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):
        """
        Initialize a :class:`.MetalComplex`.

        Parameters
        ----------
        metals : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices in :attr:`_metal_vertex_prototypes` it should
            be placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed at all :attr:`_metal_vertex_prototypes`
            on the topology graph.

        ligands : :class:`dict` or :class:`.BuildingBlock` or \
                :class:`tuple`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices in :attr:`_ligand_vertex_prototypes` it should be
            placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed at all :attr:`_ligand_vertex_prototypes`
            on the topology graph.

            If each :class:`.BuildingBlock` has a different number
            of functional groups, they can be supplied together in
            a :class:`tuple`.

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if isinstance(metals, dict):
            building_block_vertices = {
                metal: tuple(self._get_metal_vertices(ids))
                for metal, ids in metals.items()
            }
        else:
            building_block_vertices = {
                metals: tuple(self._get_metal_vertices(ids))
                for ids in range(len(self._metal_vertex_prototypes))
            }

        if isinstance(metals, dict):
            building_block_vertices.update(
                (ligand, tuple(self._get_ligand_vertices(ids)))
                for ligand, ids in ligands.items()
            )
        else:
            building_block_vertices.update(
                (ligands, tuple(self._get_metal_vertices(ids)))
                for ids in range(len(self._ligand_vertex_prototypes))
            )

        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=self._edge_prototypes,
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=None,
        )

    def _get_metal_vertices(self, vertex_ids):
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
            yield self._metal_vertex_prototypes[vertex_id]

    def _get_ligand_vertices(self, vertex_ids):
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
            yield self._ligand_vertex_prototypes[vertex_id]

    def _get_scale(self, building_block_vertices):
        return 1

    def __repr__(self):
        return (
            f'metal_complex.{self.__class__.__name__}'
            f'()'
        )
