"""
Rotaxane
========

#. :class:`.NRotaxane`

"""


import numpy as np
import rdkit.Chem.AllChem as rdkit

from .topology_graph import TopologyGraph, Vertex


class _AxleVertex(Vertex):
    """
    Places the axle in a :class:`NRotaxane`.

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

        `building_block` is placed such that its bonder-bonder
        direction vector is either parallel or anti-parallel to the
        polymer chain.

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

    def assign_func_groups_to_edges(self, building_block):
        return {}


class _CycleVertex(Vertex):
    """
    Places the cycles in a :class:`NRotaxane`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    def __init__(self, x, y, z, orientation):
        self._orientation = orientation
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

    def assign_func_groups_to_edges(self, building_block):
        return {}

    def __str__(self):
        x, y, z = self._position
        return (
            f'Vertex(id={self.id}, '
            f'position={[x, y, z]}, '
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
        repeating_unit : :class:`str`
            A string specifying the repeating unit of the macrocycles.
            For example, ``'AB'`` or ``'ABB'``. Letters are assigned to
            building block molecules in the order they are passed to
            :meth:`.ConstructedMolecule.__init__`.

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

        self._repeating_unit = repeating_unit
        self._orientations = orientations
        self._num_repeating_units = num_repeating_units

        vertices = [_AxleVertex(0, 0, 0)]
        threads = orientations * num_repeating_units
        distance = 1 / (len(threads)+1)
        for i, orientation in enumerate(threads, 1):
            vertices.append(
                _CycleVertex(
                    x=i*distance-0.5,
                    y=0,
                    z=0,
                    orientation=orientation
                )
            )
        super().__init__(tuple(vertices), (), (), num_processes)

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
        bb_map = {
            letter: bb for letter, bb in zip(threads, cycles)
        }
        building_block_vertices = {}
        building_block_vertices[axle] = self.vertices[0:1]
        for letter, vertex in zip(threads, self.vertices[1:]):
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
