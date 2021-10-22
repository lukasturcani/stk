"""
Edge
====

"""

from __future__ import annotations

from collections import abc
import typing
import numpy as np
import stk


__all__ = (
    'Edge',
)


_T = typing.TypeVar('_T', bound='Edge')


class Edge:
    """
    Represents an edge in a topology graph.

    """

    def __init__(
        self,
        id: int,
        vertex1: stk.Vertex,
        vertex2: stk.Vertex,
        periodicity: tuple[int, int, int] = (0, 0, 0),
        position: typing.Optional[np.ndarray] = None,
    ) -> None:
        """
        Initialize an :class:`.Edge` instance.

        Parameters:

            id:
                The id of the edge.

            vertex1:
                The first vertex the edge is connected to.

            vertex2:
                The second vertex the edge is connected to.

            periodicity:
                The periodicity of the edge, when going from `vertex1`
                to `vertex2`. For example, if ``(0, 0, 0)`` the edge is
                not periodic, if ``(1, 0, -1)`` the edge is periodic
                going in the positive direction along the x axis, is
                not periodic across the y axis and is periodic in the
                negative direction along the z axis.

            position:
                The position of the edge, if ``None``, the midpoint of
                the vertices is used.

        """

        if position is None:
            self._position = (
                vertex1.get_position() + vertex2.get_position()
            ) / 2
        else:
            self._position = np.array(position, dtype=np.float64)

        self._id = id
        self._vertex1_id = vertex1.get_id()
        self._vertex2_id = vertex2.get_id()
        self._periodicity = periodicity

    def get_id(self) -> int:
        """
        Get the id of the edge.

        Returns:

            The id.

        """

        return self._id

    def get_periodicity(self) -> tuple[int, int, int]:
        """
        Get the periodicity of the edge.

        Returns:

            The periodicity of the edge. For example, if ``(0, 0, 0)``
            the edge is not periodic, if ``(1, 0, -1)`` the edge is
            periodic going in the positive direction along the x axis,
            is not periodic across the y axis and is periodic in the
            negative direction along the z axis.

        """

        return self._periodicity

    def is_periodic(self) -> bool:
        """
        Return ``True`` if periodic.

        Returns:

            ``True`` if periodic.

        """

        return any(i != 0 for i in self._periodicity)

    def _with_scale(
        self: _T,
        scale: typing.Union[float, tuple[float, float, float]],
    ) -> _T:
        """
        Modify the edge.

        """

        self._position *= scale
        return self

    def with_scale(
        self,
        scale: typing.Union[float, tuple[float, float, float]],
    ) -> Edge:
        """
        Return a clone with a scaled position.

        Parameters:

            scale:
                The value by which the position of the :class:`.Edge`
                is scaled. Can be a single number if all axes are
                scaled by the same amount or a :class:`tuple` of three
                numbers if each axis is scaled by a different value.

        Returns:

            The clone.


        """

        return self.clone()._with_scale(scale)

    def clone(self) -> Edge:
        """
        Return a clone.

        Returns:

            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._id = self._id
        clone._vertex1_id = self._vertex1_id
        clone._vertex2_id = self._vertex2_id
        clone._position = np.array(self._position)
        clone._periodicity = self._periodicity
        return clone

    def get_vertex_ids(self) -> abc.Iterator[int]:
        """
        Yield the ids of the vertices connected by the edge.

        The id of the first vertex is yielded first, followed by the id
        of the second vertex.

        Yields:

            The id of a :class:`.Vertex`.

        """

        yield self._vertex1_id
        yield self._vertex2_id

    def get_vertex1_id(self) -> int:
        """
        Get the id of the first vertex.

        Returns:

            The id of the first vertex.

        """

        return self._vertex1_id

    def get_vertex2_id(self) -> int:
        """
        Get the id of the second vertex.

        Returns:

            The id of the second vertex.

        """

        return self._vertex2_id

    def get_position(self) -> np.ndarray:
        """
        Get the position.

        Returns:

            The position of the :class:`Edge`.

        """

        return np.array(self._position)

    def _with_position(
        self: _T,
        position: np.ndarray,
    ) -> _T:
        """
        Modify the edge.

        """

        self._position = np.array(position, dtype=np.float64)
        return self

    def with_position(
        self,
        position: np.ndarray,
    ) -> Edge:
        """
        Return a clone at `position`.

        Parameters:

            position:
                The desired position of the clone.

        Returns:

            The clone.

        """

        return self.clone()._with_position(position)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return (
            f'Edge({self._id}, {self._vertex1_id}, {self._vertex2_id})'
        )
