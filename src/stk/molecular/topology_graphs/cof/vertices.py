"""
COF Vertices
============

"""

import numpy as np
from scipy.spatial.distance import euclidean

from stk.utilities import get_acute_vector

from ..topology_graph import Vertex
from ..utilities import _EdgeSorter, _FunctionalGroupSorter


class _CofVertex(Vertex):
    """
    A :class:`._CofVertex` .

    """

    def __init__(self, id, position, aligner_edge=0, cell=(0, 0, 0)):
        """
        Initialize a :class:`._CofVertex`.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`tuple` of :class:`float`
            The position of the vertex.

        aligner_edge : :class:`int`, optional
            The edge which is used to align the :class:`.BuildingBlock`
            placed on the vertex. The first :class:`.FunctionalGroup`
            is rotated such that it lies exactly on this
            :class:`.Edge`. Must be between ``0`` and the number of
            edges the vertex is connected to.

        cell : :class:`tuple` of :class:`int`, optional
            The cell of the lattice in which the vertex is found.

        """

        super().__init__(id, position)
        self._aligner_edge = aligner_edge
        self._cell = np.array(cell)

    def get_aligner_edge(self):
        """
        Return the aligner edge of the vertex.

        Returns
        -------
        :class:`int`
            The aligner edge.

        """

        return self._aligner_edge

    def get_cell(self):
        return np.array(self._cell)

    @classmethod
    def init_at_center(
        cls,
        id,
        vertices,
        aligner_edge=0,
        cell=(0, 0, 0),
    ):
        """
        Initialize a :class:`._CofVertex` in the middle of `vertices`.

        Parameters
        ----------
        id : :class:`int`
            The id of the initialized vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices at whose center this one needs to be.

        aligner_edge : :class:`int`, optional
            The edge which is used to align the :class:`.BuildingBlock`
            placed on the vertex. The first :class:`.FunctionalGroup`
            is rotated such that it lies exactly on this
            :class:`.Edge`. Must be between ``0`` and the number of
            edges the vertex is connected to.

        cell : :class:`tuple` of :class:`int`, optional
            The cell of the lattice in which the vertex is found.

        Returns
        -------
        :class:`._CofVertex`
            The new vertex.

        """

        return cls(
            id=id,
            position=(
                sum(vertex.get_position() for vertex in vertices)
                / len(vertices)
            ),
            aligner_edge=aligner_edge,
            cell=cell,
        )

    @classmethod
    def init_at_shifted_center(
        cls,
        id,
        vertices,
        cell_shifts,
        lattice_constants,
        aligner_edge=0,
        cell=(0, 0, 0),
    ):
        """
        Initialize a :class:`._CofVertex` at the center of `vertices`.

        The `vertices` are shifted according to the lattice constants
        and cell shifts.

        Parameters
        ----------
        id : :class:`int`
            The id of the initialized vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices at whose center this one needs to be.

        cell_shifts : :class:`tuple` of :class:`int`
            The number of cells shifted in the x, y and z directions.

        lattice_constants : :class:`tuple` of :class:`numpy.ndarray`
            The a, b and c lattice constants.

        aligner_edge : :class:`int`, optional
            The edge which is used to align the :class:`.BuildingBlock`
            placed on the vertex. The first :class:`.FunctionalGroup`
            is rotated such that it lies exactly on this
            :class:`.Edge`. Must be between ``0`` and the number of
            edges the vertex is connected to.

        cell : :class:`tuple` of :class:`int`, optional
            The cell of the lattice in which the vertex is found.

        Returns
        -------
        :class:`._CofVertex`
            The new vertex.

        """

        positions = []
        for vertex, cell_shift in zip(vertices, cell_shifts):
            shift = sum(
                dim_shift*constant
                for dim_shift, constant
                in zip(cell_shift, lattice_constants)
            )
            positions.append(vertex.get_position() + shift)

        position = np.divide(
            np.sum(positions, axis=0),
            len(positions),
        )
        return cls(id, position, aligner_edge, cell)

    def clone(self):
        clone = super().clone()
        clone._aligner_edge = self._aligner_edge
        clone._cell = np.array(self._cell)
        return clone

    def __str__(self):
        return (
            f'Vertex(id={self._id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class LinearVertex(_CofVertex):
    def place_building_block(self, building_block, edges):
        assert (
            building_block.get_num_functional_groups() == 2
        ), (
            f'{building_block} needs to have exactly 2 functional '
            'groups but has '
            f'{building_block.get_num_functional_groups()}.'
        )
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )

        # Align the normal of the plane of best fit, defined by
        # all atoms in the building block, with the z axis.
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        normal = building_block.get_plane_normal()
        normal = get_acute_vector(
            reference=core_centroid - self._position,
            vector=normal,
        )
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=[0, 0, 1],
            origin=self._position,
        )

        # Rotate to place fg-fg vector along edge-edge vector.
        fg, = building_block.get_functional_groups(0)
        fg_centroid = building_block.get_centroid(fg.get_placer_ids())
        target = edges[0].get_position() - edges[1].get_position()
        target *= 1 if self._aligner_edge == 0 else -1

        building_block = building_block.with_rotation_between_vectors(
            start=fg_centroid - self._position,
            target=target,
            origin=self._position,
        )
        return building_block.get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg, = building_block.get_functional_groups(0)
        fg_position = building_block.get_centroid(fg.get_placer_ids())

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        edges = sorted(edges, key=fg_distance)
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class NonLinearVertex(_CofVertex):
    def place_building_block(self, building_block, edges):
        assert (
            building_block.get_num_functional_groups() > 2
        ), (
            f'{building_block} needs to have more than 2 functional '
            'groups but has '
            f'{building_block.get_num_functional_groups()}.'
        )
        # Sort to ensure that for two vertices, which are periodically
        # equivalent, "edges" has identical ordering. This means that
        # the aligner_edge is chosen consistently in both cases.
        edges = sorted(edges, key=lambda edge: edge.get_parent_id())

        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = get_acute_vector(
            reference=core_centroid - self._position,
            vector=normal,
        )
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=[0, 0, 1],
            origin=self._position,
        )
        fg, = building_block.get_functional_groups(0)
        fg_centroid = building_block.get_centroid(fg.get_placer_ids())
        edge_position = edges[self._aligner_edge].get_position()
        return building_block.with_rotation_to_minimize_angle(
            start=fg_centroid - self._position,
            target=edge_position - self._position,
            axis=np.array([0, 0, 1], dtype=np.float64),
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        # Sort to ensure that for two vertices, which are periodically
        # equivalent, "edges" has identical ordering. This means that
        # the aligner_edge is chosen consistently in both cases.
        edges = sorted(edges, key=lambda edge: edge.get_parent_id())

        fg_sorter = _FunctionalGroupSorter(building_block)
        edge_sorter = _EdgeSorter(
            edges=edges,
            aligner_edge=edges[self._aligner_edge],
            axis=fg_sorter.get_axis(),
        )
        return {
            fg_id: edge.get_id()
            for fg_id, edge in zip(
                fg_sorter.get_items(),
                edge_sorter.get_items(),
            )
        }


class UnaligningVertex(_CofVertex):
    """
    Just places a building block, does not align.

    """

    def place_building_block(self, building_block, edges):
        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):

        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }

    @classmethod
    def init_at_center(
        cls,
        id,
        vertices,
        aligner_edge=0,
        cell=(0, 0, 0),
    ):

        vertex = cls.__new__(cls)
        vertex._id = id
        vertex._position = (
            sum(vertex.get_position() for vertex in vertices)
            / len(vertices)
        )
        vertex._cell = np.array(cell)
        vertex._aligner_edge = aligner_edge
        return vertex

    @classmethod
    def init_at_shifted_center(
        cls,
        id,
        vertices,
        cell_shifts,
        lattice_constants,
        aligner_edge=0,
        cell=(0, 0, 0),
    ):
        new_vertex = cls.__new__(cls)
        new_vertex._id = id

        positions = []
        for vertex, cell_shift in zip(vertices, cell_shifts):
            shift = sum(
                dim_shift*constant
                for dim_shift, constant
                in zip(cell_shift, lattice_constants)
            )
            positions.append(vertex.get_position() + shift)

        new_vertex._position = np.divide(
            np.sum(positions, axis=0),
            len(positions),
        )
        new_vertex._cell = np.array(cell)
        new_vertex._aligner_edge = aligner_edge
        return new_vertex
