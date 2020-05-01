"""
Metal Centre
============

"""

from ..topology_graph import TopologyGraph
from ...reactions import GenericReactionFactory


class MetalCentre(TopologyGraph):
    """
    Represents a metal centre topology graph.

    This topology graph places single atoms in idealized positions
    for a metal geometry and then reacts the atoms.

    Examples
    --------
    *Construction*

    You can use :class:`.ConstructedMolecule` instances as the host,
    but you should turn them into a :class:`.BuildingBlock` first

    .. code-block:: python

        import stk

        # Build single atoms to place on metal centre topology.
        # Metal atom with 6 functional groups.
        atom = rdkit.MolFromSmiles('[Fe+2]')
        atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
        metal_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
        atom_0 = list(metal_atom.get_atoms())[0]
        metal_atom = metal_atom.with_functional_groups(
            (stk.SingleAtom(atom_0) for i in range(6))
        )


        # Nitrogen atom.
        atom = rdkit.MolFromSmiles('N')
        atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
        binding_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
        atom_0 = list(binding_atom.get_atoms())[0]
        binding_atom = binding_atom.with_functional_groups(
            (stk.SingleAtom(atom_0) for i in range(1))
        )

        # Build an Fe atom with octahedrally coordinated N atoms.
        metal_centre = stk.ConstructedMolecule(
            topology_graph=stk.metal_centre.Octahedral(
                building_blocks={
                    metal_atom: 0,
                    binding_atom: range(1, 7)
                }
            )
        )

    """

    def __init__(
        self,
        building_blocks,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):
        """
        Initialize a :class:`.MetalCentre`.

        Parameters
        ----------
        building_blocks : :class:`dict`
            :class:`dict` which maps the
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

        """

        building_blocks = {
            building_block: self._get_vertices(ids)
            for building_block, ids in building_blocks.items()
        }

        self._vertex_alignments = vertex_alignments = (
            dict(vertex_alignments)
            if vertex_alignments is not None
            else {}
        )
        super().__init__(
            building_block_vertices=building_blocks,
            edges=self._edge_prototypes,
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=None
        )

    def _get_vertices(self, vertex_ids):
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
            yield self._vertex_prototypes[vertex_id]

    def _get_scale(self, building_block_vertices):
        return 1

    def __repr__(self):
        return (
            f'metal_centre.{self.__class__.__name__}()'
        )
