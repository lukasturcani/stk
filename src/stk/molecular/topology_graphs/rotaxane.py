"""
Rotaxane
========

#. :class:`.NRotaxane`

"""


import numpy as np
import rdkit.Chem.AllChem as rdkit

from .topology_graph import TopologyGraph, VertexData, Vertex


class _AxleVertexData(VertexData):
    def get_vertex(self):
        return _AxleVertex(self)


class _AxleVertex(Vertex):
    def place_building_block(self, building_block, vertices, edges):
        building_block.set_centroid(self._position)
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return {}


class _CycleVertexData(VertexData):
    """
    Holds data for a :class:`._CycleVertex`.

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

    orientation : :class:`float`
        Can be any number from ``0`` to ``1``, both inclusive. It
        specifies the probability the building block placed on the
        vertex will have its orientation along the chain flipped.

    """

    def __init__(self, x, y, z, orientation):
        """
        Initialize a :Class:`._CycleVertexData` instance.

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

        """

        self.orientation = orientation
        super().__init__(x, y, z)

    def clone(self, clear_edges=False):
        clone = super().clone(clear_edges)
        clone.orientation = self.orientation
        return clone

    def get_vertex(self):
        return _CycleVertex(self)


class _CycleVertex(Vertex):
    """
    Places the cycles in a :class:`NRotaxane`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def __init__(self, data):
        self._orientation = data.orientation
        super().__init__(data)

    def clone(self, clear_edges=False):
        clone = super().clone(clear_edges)
        clone._orientation = self._orientation
        return clone

    def place_building_block(self, building_block, vertices, edges):
        rdkit_mol = building_block.to_rdkit_mol()
        macrocycle = max(rdkit.GetSymmSSSR(rdkit_mol), key=len)
        cycle_normal = building_block.get_plane_normal(macrocycle)
        building_block.set_centroid(
            position=self._position,
            atom_ids=macrocycle
        )
        p = [1-self._orientation, self._orientation]
        direction = np.random.choice([1, -1], p=p)
        building_block.apply_rotation_between_vectors(
            start=cycle_normal,
            target=[direction, 0, 0],
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return {}

    def __str__(self):
        return (
            f'Vertex(id={self.id}, '
            f'position={self._position.tolist()}, '
            f'orientation={self._orientation})'
        )


class NRotaxane(TopologyGraph):
    """
    Represents [n]rotaxane topology graphs.

    This class assumes one axle with (n-1) macrocycles threaded on it.
    The macrocycles are spaced evenly along the thread in repeating
    patterns. The threaded macrocycles can be described analagously
    to monomers in linear polymers, in terms of a repeating unit,
    except that no bonds are formed between them.

    The axle must be provided first to the `building_blocks` in
    :class:`.ConstructedMolecule.__init__`.

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

        cycle = stk.ConstructedMolecule(
            building_blocks=[
                stk.BuildingBlock('[Br]CC[Br]', ['bromine'])
            ],
            topology_graph=stk.macrocycle.Macrocycle('A', 5)
        )
        axle = stk.ConstructedMolecule(
            building_blocks=[
                stk.BuildingBlock('NCCN', ['amine']),
                stk.BuildingBlock('O=CCC=O', ['aldehyde'])
            ],
            topology_graph=stk.polymer.Linear('AB', 7)
        )
        rotaxane = stk.ConstructedMolecule(
            building_blocks=[axle, cycle],
            topology_graph=stk.rotaxane.NRotaxane('A', 3)
        )

    The repeating unit can also be specified through the indices of
    the building blocks

    .. code-block:: python

        # r1 and r2 are different ways to write the same thing.
        r1 = stk.ConstructedMolecule(
            building_blocks=[axle, cycle, cycle2, cycle3],
            topology_graph=stk.rotaxane.NRotaxane('ACB', 3)
        )
        r2 = stk.ConstructedMolecule(
            building_blocks=[axle, cycle, cycle2, cycle3],
            topology_graph=stk.rotaxane.NRotaxane((1, 3, 2), 3)
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
        Initialize a :class:`NRotaxane` instance.

        Parameters
        ----------
        repeating_unit : :class:`str` or :class:`tuple` of :class:`int`
            A string specifying the repeating unit of the macrocycles.
            For example, ``'AB'`` or ``'ABB'``. The first macrocycle
            passed to `building_blocks` is ``'A'`` and so on.

            The repeating unit can also be specified by the indices of
            `building_blocks`, for example ``'ABB'`` can be
            written as ``(1, 2, 2)``.

        num_repeating_units : :class:`int`
            The number of repeating units threaded along the axle.

        orientations : :class:`tuple` of :class:`float`, optional
            For each character in the repeating unit, a value
            between ``0`` and ``1`` (both inclusive) must be given in
            a :class:`tuple`. It indicates the probability that each
            macrocycle will have its orientation along the axle
            flipped. If ``0`` then the macrocycle is guaranteed not to
            flip. If ``1`` it is guaranteed to flip. This allows the
            user to create head-to-head or head-to-tail chains, as well
            as chain with a preference for head-to-head or head-to-tail
            if a number between ``0`` and ``1`` is chosen. If
            ``None`` then defaults to ``0`` in every case.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if orientations is None:
            orientations = tuple(
                0. for i in range(len(repeating_unit))
            )

        self._repeating_unit = self._normalize_repeating_unit(
            repeating_unit=repeating_unit
        )
        self._orientations = orientations
        self._num_repeating_units = num_repeating_units

        vertex_data = [_AxleVertexData(0, 0, 0)]
        threads = orientations * num_repeating_units
        distance = 1 / (len(threads)+1)
        for i, orientation in enumerate(threads, 1):
            vertex_data.append(
                _CycleVertexData(
                    x=i*distance-0.5,
                    y=0,
                    z=0,
                    orientation=orientation
                )
            )
        super().__init__(tuple(vertex_data), (), (), num_processes)

    @staticmethod
    def _normalize_repeating_unit(repeating_unit):
        if isinstance(repeating_unit, tuple):
            return repeating_unit
        base = ord('A')-1
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

        threads = self._repeating_unit*self._num_repeating_units
        axle, *cycles = building_blocks
        building_block_vertices = {}
        building_block_vertices[axle] = self.vertices[0:1]
        for bb_index, vertex in zip(threads, self.vertices[1:]):
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

        axle = next(iter(mol.building_block_vertices))
        return 0.8*axle.get_maximum_diameter()

    def __repr__(self):
        return (
            f'rotaxane.NRotaxane('
            f'{self._repeating_unit!r}, '
            f'{self._num_repeating_units}, '
            f'{self._orientations!r}'
            f')'
        )
