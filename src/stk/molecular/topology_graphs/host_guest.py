"""
Host Guest Complex
==================

#. :class:`.Complex`

"""


import numpy as np

from .topology_graph import TopologyGraph, VertexData, Vertex


class _HostVertexData(VertexData):
    def get_vertex(self):
        return _HostVertex(self)


class _HostVertex(Vertex):
    """
    Places the host in a :class:`.Complex`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def place_building_block(self, building_block, vertices, edges):
        building_block.set_centroid(self._position)
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return dict()


class _GuestVertexData(VertexData):
    """
    Holds data for a :class:`_GuestVertex`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. Must match the index in
        :attr:`TopologyGraph.vertices`.

    position : :class:`numpy.ndarray`
        The position of the vertex.

    edges : :class:`list` of :class:`.EdgeData`
        The edges connected to the vertex.

    cell : :class:`numpy.ndarray`
        The unit cell in which the vertex is found.

    start : :class:`numpy.ndarray`
        A direction vector which gets aligned with `target`.

    target : :class:`numpy.ndarray`
        A direction vector which determines the rotation applied to
        the placed building block. A rotation such that
        `start` is transformed into `target` is applied
        to the placed building block.

    """

    def __init__(self, x, y, z, start, target):
        """
        Initialize a :class:`_GuestVertexData` instance.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        start : :class:`numpy.ndarray`
            A direction vector which gets aligned with `target`.

        target : :class:`numpy.ndarray`
            A direction vector which determines the rotation applied to
            the placed building block. A rotation such that
            `start` is transformed into `target` is applied
            to the placed building block.

        """

        self.start = start
        self.target = target
        super().__init__(x, y, z)

    def clone(self, clear_edges=False):
        clone = super().clone(clear_edges)
        clone.start = np.array(self.start)
        clone.target = np.array(self.target)
        return clone

    def get_vertex(self):
        return _GuestVertex(self)


class _GuestVertex(Vertex):
    """
    Places the guest in a :class:`.Complex`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def __init__(self, data):
        self._start = np.array(data.start)
        self._target = np.array(data.target)
        super().__init__(data)

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            If ``True`` the :attr:`edges` attribute of the clone will
            be empty.

        Returns
        -------
        :class:`Vertex`
            The clone.

        """

        clone = super().clone(clear_edges)
        clone._start = np.array(self._start)
        clone._target = np.array(self._target)
        return clone

    def place_building_block(self, building_block, vertices, edges):
        building_block.set_centroid(self._position)
        building_block.apply_rotation_between_vectors(
            start=self._start,
            target=self._target,
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return dict()


class Complex(TopologyGraph):
    """
    Represents a host-guest complex topology graph.

    When using this topology graph, the host must be first in the
    `building_blocks` of the :class:`.ConstructedMolecule`
    and the guest must be second.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    Examples
    --------
    .. code-block:: python

        import stk

        guest = stk.BuildingBlock('[Br][Br]')
        host = stk.ConstructedMolecule(
            building_blocks=[
                stk.BuildingBlock('NCCN', ['amine']),
                stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
            ],
            topology_graph=stk.cage.FourPlusSix()
        )
        complex1 = stk.ConstructedMolecule(
            building_blocks=[host, guest],
            topology_graph=stk.host_guest.Complex()
        )

    Change the position and orientation of the guest

    .. code-block:: python

        complex2 = stk.ConstructedMolecule(
            building_blocks=[host, guest],
            topology_graph=stk.host_guest.Complex(
                # Apply a rotation onto the guest molecule such that
                # the vector returned by get_direction() has the same
                # direction as [1, 1, 1].
                guest_start=guest.get_direction(),
                guest_target=[1, 1, 1],
                displacement=[5.3, 2.1, 7.1]
            )
        )

    """

    def __init__(
        self,
        guest_start=None,
        guest_target=None,
        displacement=None,
        num_processes=1
    ):
        """
        Initialize an instance of :class:`.Complex`.

        Parameters
        ----------
        guest_start : :class:`tuple` of :class:`float`, optional
            A direction vector which gets aligned with `guest_target`.

        guest_target : :class:`tuple` of :class:`float`, optional
            A direction vector which determines the rotation applied to
            the guest building block. A rotation such that
            `guest_start` is transformed into `guest_target` is applied
            to the guest building block.

        displacement : :class:`list` of :class:`float`, optional
            The translational offset of the guest from the center of
            the host cavity.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        Raises
        ------
        :class:`TypeError`
            If `guest_start` or `guest_target` is defined but the other
            is not.

        """

        num_nones = sum(
            1 for vector in (guest_start, guest_target)
            if vector is None
        )
        if num_nones == 1:
            raise TypeError(
                'If guest_start or guest_target is defined, '
                'the other must be too.'
            )

        if guest_start is not None:
            start = guest_start = tuple(guest_start)
            target = guest_target = tuple(guest_target)
        else:
            start = target = (1., 0., 0.)

        # Save the values as None, for __repr__.
        self._guest_start = guest_start
        self._guest_target = guest_target
        self._displacement = displacement

        if displacement is None:
            displacement = np.array([0, 0, 0])

        x, y, z = displacement
        vertices = (
            _HostVertexData(0, 0, 0),
            _GuestVertexData(x, y, z, start, target)
        )
        super().__init__(vertices, (), (), num_processes)

    def assign_building_blocks_to_vertices(self, building_blocks):
        """
        Assign `building_blocks` to :attr:`vertices`.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        Returns
        -------
        :class:`dict`
            Maps the `building_blocks`, to the
            :class:`~.topologies.base.Vertex` objects in
            :attr:`vertices` they are placed on during construction.
            The :class:`dict` has the form

            .. code-block:: python

                building_block_vertices = {
                    BuildingBlock(...): [Vertex(...), Vertex(...)],
                    BuildingBlock(...): [
                        Vertex(...),
                        Vertex(...),
                        Vertex(...),
                    ]
                    ConstructedMolecule(...): [Vertex(...)]
                }

        """

        return {
            bb: [vertex]
            for bb, vertex in zip(building_blocks, self.vertices)
        }

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

        return 1

    def __repr__(self):
        return (
            f'host_guest.Complex('
            f'guest_start={self._guest_start!r}, '
            f'guest_target={self._guest_target!r}, '
            f'displacement={self._displacement!r})'
        )
