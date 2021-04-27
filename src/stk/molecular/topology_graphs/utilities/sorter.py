"""
Sorter
======

"""

import numpy as np

from stk.utilities import vector_angle


class _Sorter:
    """
    Sorts items according to their angle from a reference vector.

    """

    __slots__ = ['_items', '_reference', '_axis']

    def __init__(self, items, reference, axis):
        """
        Initialize a :class:`._Sorter`.

        Parameters
        ----------
        items : :class:`iterable` of :class:`object`
            The items to be sorted.

        reference : :class:`numpy.ndarray`
            The reference from which the angle is calculated.

        axis : :class:`numpy.ndarray`
            A vector orthogonal to `reference`, used to determine
            which direction is clockwise. Must be an immutable
            array.

        """

        self._items = items
        self._reference = reference
        self._axis = axis

    def _get_vector(self, item):
        """
        Get the vector according to which `item` should be sorted.

        Parameters
        ----------
        item : :class:`object`
            The item being sorted.

        Returns
        -------
        :class:`numpy.ndarray`
            The vector of the `item`, which should be used to measure
            its angle with respect to the *reference*.

        """

        raise NotImplementedError()

    def _get_angle(self, item):
        """
        Get the angle of `vector` relative to `reference`.

        Parameters
        ----------
        item : :class:`object`
            The item being sorted.

        Returns
        -------
        :class:`float`
            The angle between `item` and the reference vector.

        """

        vector = self._get_vector(item)
        theta = vector_angle(self._reference, vector)
        projection = vector @ self._axis
        if theta > 0 and projection < 0:
            return 2*np.pi - theta
        return theta

    def get_items(self):
        """
        Yield the sorted items.

        Yields
        ------
        :class:`object`
            An item.

        """

        yield from sorted(self._items, key=self._get_angle)

    def get_axis(self):
        """
        Get the axis used to determine which direction is clockwise.

        Returns
        -------
        :class:`numpy.ndarray`
            The axis. The array is immutable.

        """

        return self._axis
