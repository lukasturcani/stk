import numpy as np
import rdkit.Chem.AllChem as rdkit
from collections import defaultdict

from ..topology_graph import TopologyGraph, Vertex


class CageVertex(Vertex):
    """
    Represents a vertex of a :class:`.CageTopology`.

    Attributes
    ----------
    aligner_edge : :class:`.Edge`
        The :class:`.Edge` in :attr:`edges`, which is used to align the
        :class:`.BuildingBlock` placed on the vertex. The first
        :class:`.FunctionalGroup` in :attr:`.BuildingBlock.func_groups`
        is rotated such that it lies exactly on this :class:`.Edge`.

    Methods
    -------
    :meth:`__init__`
    :meth:`clone`
    :meth:`place_building_block`

    """

    def __init__(self, x, y, z):
        self.aligner_edge = self.edges[0]
        super().__init__(x, y, z)

    def clone(self, clear_edges=False):
        """
        Create a clone of the instance.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            If ``True`` the :attr:`edges` attribute of the clone will
            be empty.

        Returns
        -------
        :class:`Vertex`
            A clone with the same position but not connected to any
            :class:`.Edge` objects.

        """

        clone = super().clone(clear_edges)
        clone.alignment = self.alignment
        return clone

    def place_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        if len(building_block.func_groups) == 2:
            return self._place_linear_building_block(building_block)
        return self._place_nonlinear_building_block(building_block)

    def _place_linear_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(
            position=self._coord,
            atom_ids=building_block.get_bonder_ids()
        )
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_atom_ids()
        )
        start = fg_centroid - self._coord
        edge_coord = self.aligner_edge.get_position()
        target = edge_coord - self._get_edge_centroid()
        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=self._coord
        )

    def _place_nonlinear_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(
            position=self._coord,
            atom_ids=building_block.get_bonder_ids()
        )
        edge_normal = self._get_edge_plane_normal()
        building_block.apply_rotation_between_vectors(
            start=building_block.get_bonder_plane_normal(),
            target=edge_normal,
            origin=self._coord
        )
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_atom_ids()
        )
        start = fg_centroid - self._coord
        edge_coord = self.aligner_edge.get_position()
        target = edge_coord - self._get_edge_centroid()
        building_block.apply_rotation_to_minimize_theta(
            start=start,
            target=target,
            axis=edge_normal,
            origin=self._coord
        )
        return building_block.get_position_matrix()

    def _assign_func_groups_to_edges(self, building_block, fg_map):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be associated with one of the :class:`.Edge` instances in
        :attr:`edges`. Then, using `fg_map`, the
        :class:`FunctionalGroup` instances in the molecule being
        constructed need to be assigned to those edges. This is
        because bonds need to be formed between functional groups of
        the molecule being constructed, not the `building_block`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to

        fg_map : :class:`dict`
            A mapping from :class:`.FunctionalGroup` instances in
            `building_block` to the equivalent
            :class:`.FunctionalGroup` instances in the molecule being
            constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        ...


class CageTopology(TopologyGraph):
    """
    A base :class:`.TopologyGraph` for cages.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.CageVertex`
        A class attribute. It holds vertices used to make a specific
        cage topology graph. This needs to be defined by a subclass.

    edges : :class:`tuple` of :class:`.Edge`
        A class attribute. It hold the edges used to make a specific
        cage topology graph. This needs to be defined by a subclass.

    """

    def __init__(self, vertex_alignments=None, processes=1):
        """
        Initialize a :class:`.CageTopology`.

        Parmaeters
        ----------
        vertex_alignments : :class:`dict`, optional
            A mapping from a :class:`.CageVertex` in :attr:`vertices`
            to an :class:`.Edge` connected to it. The :class:`Edge` is
            used to align the first :class:`.FunctionalGroup` of a
            :class:`.BuildingBlock` placed on that vertex. Only
            vertices which need to have their default edge changed need
            to be present in the :class:`dict`. If ``None`` then the
            first :class:`.Edge` in :class:`.Vertex.edges` is for each
            vertex is used. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers.

        processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if vertex_alignments is None:
            vertex_alignments = {}

        vertices = []
        for vertex in self.vertices:
            clone = vertex.clone()
            edge = vertex_alignments.get(vertex, 0)
            clone.aligner_edge = edge
            vertices.append(clone)
        super().__init__(tuple(vertices), self.edges, processes)

    def _assign_building_blocks_to_vertices(
        self,
        mol,
        building_blocks
    ):
        """
        Assign `building_blocks` to :attr:`vertices`.

        This method will assign a random building block with the
        correct amount of functional groups to each vertex.

        Note
        ----
        This method will modify
        :attr:`.ConstructedMolecule.building_block_vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance being
            constructed.

        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        Returns
        -------
        None : :class:`NoneType`

        """

        bb_by_degree = defaultdict(list)
        for bb in building_blocks:
            bb_by_degree[len(bb.func_groups)].append(bb)

        for vertex in self.vertices:
            bb = np.random.choice(bb_by_degree[len(vertex.edges)])
            mol.building_block_vertices[bb].append(vertex)

    def _get_scale(self, mol):
        """
        Get the scale used for the positions of :attr:`vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        :class:`float` or :class:`list` of :class:`float`
            The value by which the position of each :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        """

        return max(
            bb.get_maximum_diameter()
            for bb in mol.building_block_vertices
        )
