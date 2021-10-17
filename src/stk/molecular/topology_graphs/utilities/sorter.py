"""
Sorter
======

"""

import typing
import numpy as np

from collections import abc
from stk.utilities import vector_angle


__all__ = (
    'Sorter',
    'IHasPosition',
)


class IHasPosition(typing.Protocol):
    def get_position(self) -> np.ndarray:
        pass


_T = typing.TypeVar('_T', bound=IHasPosition)


class Sorter(typing.Generic[_T]):
    """
    Sorts items according to their angle from a reference vector.

    """

    __slots__ = ['_items', '_reference', '_axis']

    def __init__(
        self,
        items: abc.Iterable[_T],
        reference: np.ndarray,
        axis: np.ndarray,
    ) -> None:
        """
        Initialize a :class:`.Sorter`.

        Parameters:

            items:
                The items to be sorted.

            reference:
                The reference from which the angle is calculated.

            axis:
                A vector orthogonal to `reference`, used to determine
                which direction is clockwise. Must be an immutable
                array.

        """

        self._items = items
        self._reference = reference
        self._axis = axis

    def _get_vector(
        self,
        item: _T,
    ) -> np.ndarray:
        """
        Get the vector according to which `item` should be sorted.

        Parameters:

            item:
                The item being sorted.

        Returns:

            The vector of the `item`, which should be used to measure
            its angle with respect to the *reference*.

        """

        raise NotImplementedError()

    def _get_angle(
        self,
        item: _T,
    ) -> float:
        """
        Get the angle of `vector` relative to `reference`.

        Parameters:

            item: The item being sorted.

        Returns:

            The angle between `item` and the reference vector.

        """

        vector = self._get_vector(item)
        theta = vector_angle(self._reference, vector)
        projection = vector @ self._axis
        if theta > 0 and projection < 0:
            return 2*np.pi - theta
        return theta

    def get_items(self) -> abc.Iterator[_T]:
        """
        Yield the sorted items.

        Yields:

            An item.

        """

        yield from sorted(self._items, key=self._get_angle)

    def get_axis(self) -> np.ndarray:
        """
        Get the axis used to determine which direction is clockwise.

        Returns:

            The axis. The array is immutable.

        """

        return self._axis
