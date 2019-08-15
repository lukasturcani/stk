import numpy as np

from .topology_graph import TopologyGraph, Vertex


class _HostVertex(Vertex):
    """
    Places the host in a :class:`.Complex`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    def place_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(self._position)
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(self, building_block, fg_map):
        return


class _GuestVertex(Vertex):
    """
    Places the guest in a :class:`.Complex`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    def __init__(self, x, y, z, axis, angle):
        self._axis = axis
        self._angle = angle
        super().__init__(x, y, z)

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
        clone._axis = self._axis
        clone._angle = self._angle
        return clone

    def place_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(self._position)
        if self._axis is not None:
            building_block.apply_rotation_about_axis(
                angle=self._angle,
                axis=self._axis,
                origin=self._position
            )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(self, building_block, fg_map):
        return


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
            topology_graph=stk.host_guest_complex.Complex()
        )

    Change the position and orientation of the guest

    .. code-block:: python

        complex2 = stk.ConstructedMolecule(
            building_blocks=[host, guest],
            topology_graph=stk.host_guest_complex.Complex(
                axis=[1, 0, 0],
                angle=3.14/3,
                displacement=[5.3, 2.1, 7.1]
            )
        )


    """

    def __init__(
        self,
        axis=None,
        angle=0,
        displacement=None,
        processes=1
    ):
        """
        Initialize an instance of :class:`.Complex`.

        Parameters
        ----------
        axis : :class:`list` of :class:`int`, optional
            The axis about which the guest is rotated.

        angle : :class:`float`, optional
            The angle by which the guest is rotated in radians.

        displacement : :class:`list` of :class:`float`, optional
            The translational offset of the guest from the center of
            the host cavity.

        processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        self._axis = axis
        self._angle = angle
        self._displacement = displacement

        if displacement is None:
            displacement = np.array([0, 0, 0])

        x, y, z = displacement
        vertices = (
            _HostVertex(0, 0, 0),
            _GuestVertex(x, y, z, axis, angle)
        )
        super().__init__(vertices, (), processes)

    def _assign_building_blocks_to_vertices(
        self,
        mol,
        building_blocks
    ):
        """
        Assign `building_blocks` to :attr:`vertices`.

        Assignment is done by modifying
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

        host, guest = building_blocks
        mol.building_block_vertices[host].append(self.vertices[0])
        mol.building_block_vertices[guest].append(self.vertices[1])

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
            f'host_guest_complex.Complex(axis={self._axis!r}, '
            f'angle={self._angle!r}, '
            f'displacement={self._displacement!r})'
        )
