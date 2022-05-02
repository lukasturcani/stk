"""
Vertex
======

"""

from __future__ import annotations

import typing

import numpy as np

from ...molecules import BuildingBlock
from .edge import Edge

_VertexT = typing.TypeVar('_VertexT', bound='Vertex')


class Vertex:
    """
    An abstract base class for :class:`.TopologyGraph` vertices.

    Notes:

        You might notice that some of the public methods of this
        abstract base class are implemented. This is purely for
        convenience when implementing subclasses. The implemented
        public methods are simply default implementations, which can
        be safely ignored or overridden, when implementing subclasses.
        Any private methods are implementation details of these
        default implementations.

    """

    def __init__(
        self,
        id: int,
        position: typing.Union[tuple[float, float, float], np.ndarray],
    ) -> None:
        """
        Initialize a :class:`.Vertex`.

        Parameters:

            id:
                The id of the vertex.

            position:
                The position of the vertex.

        """

        self._id = id
        self._position = np.array(position, dtype=np.float64)

    def get_id(self) -> int:
        """
        Get the id.

        Returns:

            The id.

        """

        return self._id

    def _with_scale(
        self: _VertexT,
        scale: typing.Union[float, tuple[float, float, float]],
    ) -> _VertexT:
        """
        Modify the vertex.

        """

        self._position *= scale
        return self

    def with_scale(
        self,
        scale: typing.Union[float, tuple[float, float, float]],
    ) -> Vertex:
        """
        Get a clone with a scaled position.

        Parameters:

            scale:
                The value by which the position of the :class:`Vertex`
                is scaled. Can be a single number if all axes are
                scaled by the same amount or a :class:`tuple` of three
                numbers if each axis is scaled by a different value.

        Returns:

            The clone.

        """

        return self.clone()._with_scale(scale)

    def clone(self) -> Vertex:
        """
        Return a clone.

        Returns:

            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        Vertex.__init__(clone, self._id, self._position)
        return clone

    def get_position(self) -> np.ndarray:
        """
        Get the position.

        Returns:

            The position of the :class:`Vertex`.

        """

        return np.array(self._position)

    def _with_position(
        self: _VertexT,
        position: np.ndarray,
    ) -> _VertexT:
        """
        Modify the vertex.

        """

        self._position = np.array(position, dtype=np.float64)
        return self

    def with_position(
        self,
        position: np.ndarray,
    ) -> Vertex:
        """
        Get a clone at a certain position.

        Parameters:

            position:
                The desired position of the clone.

        Returns:

            The clone.

        """

        return self.clone()._with_position(position)

    def get_cell(self) -> np.ndarray:
        """
        Get the cell of the lattice in which the vertex is found.

        Returns:

            The cell of the lattice in which the vertex is found.

        """

        return np.array([0, 0, 0])

    def place_building_block(
        self,
        building_block: BuildingBlock,
        edges: tuple[Edge, ...],
    ) -> np.ndarray:
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters:

            building_block:
                The building block molecule which is to be placed on
                the vertex.

            edges:
                The edges to which the vertex is attached.

        Returns:

            The position matrix of `building_block` after being
            placed.

        """

        raise NotImplementedError()

    def map_functional_groups_to_edges(
        self,
        building_block: BuildingBlock,
        edges: tuple[Edge, ...],
    ) -> dict[int, int]:
        """
        Map functional groups to edges.

        Each functional group in `building_block` needs to be assigned
        to an edge in `edges`.

        Parameters:

            building_block : :class:`.BuildingBlock`
                The building block which is needs to have functional
                groups assigned to edges.

            edges:
                The edges to which the vertex is attached.

        Returns:

            A mapping from the id of a functional group in
            `building_block` to the id of the edge in `edges` it
            is assigned to.

        """

        raise NotImplementedError()

    def __str__(self) -> str:
        position = self._position.tolist()
        return f'Vertex(id={self._id}, position={position})'

    def __repr__(self) -> str:
        return str(self)
