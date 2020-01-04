import numpy as np
from .edge import Edge


class Edge_(Edge):
    """
    An implementation of the :class:`.Edge` interface.

    """

    def __init__(self, data):
        """
        Initialize an :class:`.Edge_`.

        Parameters
        ----------
        data : :class:`.EdgeData`
            The edge data.

        """

        self._vertex_ids = tuple(
            vertex.get_id() for vertex in data.get_vertices()
        )
        self._id = data.get_id()
        self._periodicity = data.get_periodicity()

        self._has_custom_position = data.has_custom_position()
        self._position = data.get_position()
        self._lattice_constants = tuple(data.get_lattice_constants())

    def get_periodicity(self):
        return np.array(self._periodicity)

    def _with_scale(self, scale):
        """
        Modify the edge.

        """

        self._position *= scale
        self._lattice_constants = tuple(
            scale*constant for constant in self._lattice_constants
        )
        return self

    def with_scale(self, scale):
        return self.clone()._with_scale(scale)

    def clone(self):
        clone = self.__class__.__new__(self.__class__)
        clone._id = self._id
        clone._has_custom_position = self._has_custom_position
        clone._periodicity = np.array(self._periodicity)
        clone._lattice_constants = tuple(
            np.array(constant) for constant in self._lattice_constants
        )
        clone._vertex_ids = self._vertex_ids
        clone._position = np.array(self._position)
        return clone

    def get_vertex_ids(self):
        yield from self._vertex_ids

    def get_position(self, reference=None, vertices=None):
        if reference is None or not self.is_periodic():
            return np.array(self._position)

        other = vertices[
            next(v for v in self._vertex_ids if v != reference.id)
        ]
        direction = (
            1 if reference is vertices[self._vertex_ids[0]] else -1
        )
        end_cell = reference.get_cell() + direction*self._periodicity
        cell_shift = end_cell - other.get_cell()
        shift = 0
        for dim, constant in zip(cell_shift, self._lattice_constants):
            shift += dim*constant
        return (
            (other.get_position()+shift+reference.get_position()) / 2
        )

    def __str__(self):
        return repr(self)

    def __repr__(self):
        vertices = ', '.join(str(id_) for id_ in self._vertex_ids)
        if self._custom_position:
            position = f', position={self._position!r}'
        else:
            position = ''

        if any(i != 0 for i in self._periodicity):
            periodicity = f', periodicity={tuple(self._periodicity)!r}'
        else:
            periodicity = ''

        return f'Edge({vertices}{position}{periodicity})'
