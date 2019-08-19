"""
Macrocycle
==========

#. :class:`.Macrocycle`

"""

import numpy as np
import logging
from scipy.spatial.distance import euclidean

from .topology_graph import TopologyGraph, Vertex, Edge


logger = logging.getLogger(__name__)


class _CycleVertex(Vertex):
    """
    Represents a vertex in the middle of a linear polymer chain.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    def __init__(self, x, y, z, orientation, angle):
        """
        Initialize a :class:`.LinearVertex`.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        orientation : :class:`float`
            Can be any number from ``0`` to ``1``, both inclusive. It
            specifies the probability the building block placed on the
            vertex will have its orientation along the chain flipped.

        angle : :class:`float`
            The angle along the macrocycle at which the vertex is
            found.

        """

        self._orientation = orientation
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
        clone._orientation = self._orientation
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

        if len(building_block.func_groups) > 2:
            logger.warning(
                'You are placing a building block which has more than '
                'two functional groups along the backbone of '
                'a Macrocycle topology. You can remove extra '
                'functional groups from the func_groups attribute to '
                'remove this message.'
            )

        building_block.set_centroid(
            position=self._position,
            atom_ids=building_block.get_bonder_ids(fg_ids=(0, 1))
        )
        bonder_vector = next(
            building_block.get_bonder_direction_vectors(
                fg_ids=(0, 1)
            )
        )[-1]

        p = [1-self._orientation, self._orientation]
        direction = np.random.choice([1, -1], p=p)
        building_block.apply_rotation_between_vectors(
            start=bonder_vector,
            target=[direction, 0, 0],
            origin=self._position
        )
        building_block.apply_rotation_about_axis(
            angle=self._angle-(np.pi/2),
            axis=np.array([0, 0, 1]),
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(self, building_block):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be associated with one of the :class:`.Edge` instances in
        :attr:`edges`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        Returns
        -------
        :class:`dict`
            A mapping from the id of a functional group in
            `building_block` to the id of the edge in :attr:`edges` it
            is assigned to.

        """

        return {
            fg_id: e.id for fg_id, e in enumerate(sorted(
                self.edges,
                key=lambda e: self._fg0_distance(building_block, e)
            ))
        }

    def _fg0_distance(self, building_block, e):
        fg_position = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        return euclidean(e.get_position(), fg_position)

    def __str__(self):
        x, y, z = self._position
        return (
            f'Vertex(id={self.id}, '
            f'position={[x, y, z]}, '
            f'orientation={self._orientation}, '
            f'angle={self._angle})'
        )


class Macrocycle(TopologyGraph):
    """
    Represents macrocycle topology graphs.

    The macrocycle can be represented as a linear polymer with the two
    end groups bonded to close the loop.

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

        macrocycle = stk.ConstructedMolecule(
            building_blocks=[
                stk.BuildingBlock('NCCN', ['amine']),
                stk.BuildingBlock('O=CCC=O', ['aldehyde'])
            ],
            topology_graph=stk.macrocycle.Macrocycle('AB', 5)
        )

    """

    def __init__(
        self,
        repeating_unit,
        num_repeating_units,
        orientations=None,
        num_processes=1
    ):
        """
        Initialize a :class:`Macrocycle` instance.

        Parameters
        ----------
        repeating_unit : :class:`str`
            A string specifying the repeating unit of the macrocycle.
            For example, ``'AB'`` or ``'ABB'``. Letters are assigned to
            building block molecules in the order they are passed to
            :meth:`.ConstructedMolecule.__init__`.

        num_repeating_units : :class:`int`
            The number of repeating units which are used to make the
            macrocycle.

        orientations : :class:`tuple` of :class:`float`, optional
            For each character in the repeating unit, a value
            between ``0`` and ``1`` (both inclusive) must be given in
            a :class:`tuple`. It indicates the probability that each
            monomer will have its orientation along the chain flipped.
            If ``0`` then the monomer is guaranteed not to flip. If
            ``1`` it is guaranteed to flip. This allows the user to
            create head-to-head or head-to-tail chains, as well as
            chain with a preference for head-to-head or head-to-tail if
            a number between ``0`` and ``1`` is chosen. If ``None``
            then ``0`` is picked in every case.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if orientations is None:
            orientations = tuple(
                0. for i in range(len(repeating_unit))
            )

        # Keep these for __repr__
        self._repeating_unit = repeating_unit
        self._orientations = orientations
        self._num_repeating_units = num_repeating_units

        chain = orientations*num_repeating_units
        # Each monomer in the macrocycle is separated by angle_diff.
        angle_diff = (2*np.pi)/len(chain)
        vertices = []
        edges = []
        for i, orientation in enumerate(chain):
            theta = i*angle_diff
            v = _CycleVertex(
                x=np.cos(theta),
                y=np.sin(theta),
                z=0,
                orientation=orientation,
                angle=theta
            )
            vertices.append(v)

            if i > 0:
                edges.append(Edge(vertices[i-1], vertices[i]))

        edges.append(Edge(vertices[0], vertices[-1]))
        super().__init__(
            vertices=tuple(vertices),
            edges=tuple(edges),
            construction_stages=(),
            num_processes=num_processes
        )

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

        polymer = self._repeating_unit*self._num_repeating_units
        bb_map = {
            letter: bb for letter, bb in zip(polymer, building_blocks)
        }
        building_block_vertices = {}
        for letter, vertex in zip(polymer, self.vertices):
            bb = bb_map[letter]
            building_block_vertices[bb] = (
                building_block_vertices.get(bb, [])
            )
            building_block_vertices[bb].append(vertex)
        return building_block_vertices

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

        length = len(self._repeating_unit)*self._num_repeating_units
        return length*0.25*max(
            bb.get_maximum_diameter()
            for bb in mol.building_block_vertices
        )

    def __repr__(self):
        return (
            f'macrocycle.Macrocycle({self._repeating_unit!r}, '
            f'{self._num_repeating_units!r}, '
            f'{self._orientations!r})'
        )
