"""
Macrocycle
==========

#. :class:`.Macrocycle`

"""

import numpy as np
import logging
from scipy.spatial.distance import euclidean

from .topology_graph import TopologyGraph, VertexData, Vertex, EdgeData


logger = logging.getLogger(__name__)


class _CycleVertexData(VertexData):
    """
    Holds data for a :class:`.CycleVertex`.

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

    flip : :class:`bool`
        If ``True`` any building block placed by the vertex will
        have its orientation along the chain flipped.

    angle : :class:`float`
        The angle along the macrocycle at which the vertex is
        found.

    """

    def __init__(self, x, y, z, flip, angle):
        """
        Initialize a :class:`._CycleVertexData`.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

    flip : :class:`bool`
        If ``True`` any building block placed by the vertex will
        have its orientation along the chain flipped.

        angle : :class:`float`
            The angle along the macrocycle at which the vertex is
            found.

        """

        self.flip = flip
        self.angle = angle
        super().__init__(x, y, z)

    def clone(self, clear_edges=False):
        clone = super().clone(clear_edges)
        clone.flip = self.flip
        clone.angle = self.angle
        return clone

    def get_vertex(self):
        return _CycleVertex(self)


class _CycleVertex(Vertex):
    """
    Represents a vertex in the middle of a linear polymer chain.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def __init__(self, data):
        self._flip = data.flip
        self._angle = data.angle
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
        clone._flip = self._flip
        clone._angle = self._angle
        return clone

    def place_building_block(self, building_block, vertices, edges):
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
        building_block.apply_rotation_between_vectors(
            start=bonder_vector,
            target=[-1 if self._flip else 1, 0, 0],
            origin=self._position
        )
        building_block.apply_rotation_about_axis(
            angle=self._angle-(np.pi/2),
            axis=np.array([0, 0, 1]),
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return {
            fg_id: edge_id for fg_id, edge_id in enumerate(sorted(
                self._edge_ids,
                key=lambda edge_id: self._fg0_distance(
                    building_block=building_block,
                    edge_id=edge_id,
                    edges=edges
                )
            ))
        }

    def _fg0_distance(self, building_block, edge_id, edges):
        fg_position = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        return euclidean(edges[edge_id].get_position(), fg_position)

    def __str__(self):
        return (
            f'Vertex(id={self.id}, '
            f'position={self._position.tolist()}, '
            f'flip={self._flip}, '
            f'angle={self._angle})'
        )


class Macrocycle(TopologyGraph):
    """
    Represents a macrocycle topology graph.

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

    The repeating unit can also be specified through the indices of
    the building blocks

    .. code-block:: python

        bb1 = stk.BuildingBlock('BrCCBr', ['bromine'])
        bb2 = stk.BuildingBlock('BrCNCBr', ['bromine'])
        bb3 = stk.BuildingBlock('BrCNNCBr', ['bromine'])

        # c1 and c2 are different ways to write the same thing.
        c1 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2, bb3],
            topology_graph=stk.macrocycle.Macrocycle('ACB', 3)
        )
        c2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2, bb3],
            topology_graph=stk.macrocycle.Macrocycle((0, 2, 1), 3)
        )

    :class:`.Macrocycle` shares many parameters with :class:`.Linear`,
    and the examples described there are also valid for this class.
    Be sure to read them.

    """

    def __init__(
        self,
        repeating_unit,
        num_repeating_units,
        orientations=None,
        random_seed=None,
        num_processes=1
    ):
        """
        Initialize a :class:`Macrocycle` instance.

        Parameters
        ----------
        repeating_unit : :class:`str` or :class:`tuple` of :class:`int`
            A string specifying the repeating unit of the macrocycle.
            For example, ``'AB'`` or ``'ABB'``. The first building
            block passed to `building_blocks` is ``'A'`` and so on.

            The repeating unit can also be specified by the indices of
            `building_blocks`, for example ``'ABB'`` can be
            written as ``(0, 1, 1)``.

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

            It is also possible to supply an orientation for every
            vertex in the final topology graph. In this case, the
            length of `orientations` must be equal to
            ``len(repeating_unit)*num_repeating_units``.

        random_seed : :class:`int`, optional
            The random seed to use when choosing random orientations.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if orientations is None:
            orientations = tuple(
                0. for i in range(len(repeating_unit))
            )

        if len(orientations) == len(repeating_unit):
            orientations = orientations*num_repeating_units

        chain_length = len(repeating_unit)*num_repeating_units
        if len(orientations) != chain_length:
            raise ValueError(
                'The length of orientations must match either '
                'the length of repeating_unit or the '
                'total number of vertices.'
            )

        generator = np.random.RandomState(random_seed)

        # Keep these for __repr__.
        self._repeating_unit = self._normalize_repeating_unit(
            repeating_unit=repeating_unit
        )
        self._num_repeating_units = num_repeating_units

        # Each monomer in the macrocycle is separated by angle_diff.
        angle_diff = (2*np.pi)/chain_length
        vertex_data = []
        edge_data = []
        choices = [True, False]
        for i, p in enumerate(orientations):
            theta = i*angle_diff
            v = _CycleVertexData(
                x=np.cos(theta),
                y=np.sin(theta),
                z=0,
                flip=generator.choice(choices, p=[p, 1-p]),
                angle=theta
            )
            vertex_data.append(v)

            if i > 0:
                edge_data.append(
                    EdgeData(vertex_data[i-1], vertex_data[i])
                )

        # Save the chosen orientations for __repr__.
        self._orientations = tuple(int(v.flip) for v in vertex_data)

        edge_data.append(EdgeData(vertex_data[0], vertex_data[-1]))
        super().__init__(
            vertex_data=tuple(vertex_data),
            edge_data=tuple(edge_data),
            construction_stages=(),
            num_processes=num_processes
        )

    @staticmethod
    def _normalize_repeating_unit(repeating_unit):
        if isinstance(repeating_unit, tuple):
            return repeating_unit
        base = ord('A')
        return tuple(ord(letter)-base for letter in repeating_unit)

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
        building_block_vertices = {}
        for bb_index, vertex in zip(polymer, self.vertices):
            bb = building_blocks[bb_index]
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
