"""
Macrocycle
==========

#. :class:`.Macrocycle`

"""

import numpy as np

from .vertices import _CycleVertex
from ..topology_graph import TopologyGraph, Edge
from ...reactions import GenericReactionFactory


class Macrocycle(TopologyGraph):
    """
    Represents a macrocycle topology graph.

    Examples
    --------
    *Construction*

    This topology graph essentially makes a polymer chain and joins
    the ends, hence the constructor parameters allows you to specify
    the chain

    .. code-block:: python

        import stk

        macrocycle = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='AB',
                num_repeating_units=5,
            ),
        )

    The repeating unit can also be specified through the indices of
    the building blocks

    .. code-block:: python

        bb1 = stk.BuildingBlock('BrCCBr', ['bromine'])
        bb2 = stk.BuildingBlock('BrCNCBr', ['bromine'])
        bb3 = stk.BuildingBlock('BrCNNCBr', ['bromine'])

        # c1 and c2 are different ways to write the same thing.
        c1 = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(bb1, bb2, bb3),
                repeating_unit='ACB',
                num_repeating_units=3
            )
        )
        c2 = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(bb1, bb2, bb3),
                repeating_unit=(0, 2, 1),
                num_repeating_units=3,
            )
        )

    """

    def __init__(
        self,
        building_blocks,
        repeating_unit,
        num_repeating_units,
        orientations=None,
        random_seed=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1
    ):
        """
        Initialize a :class:`Macrocycle` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks of the macrocycle.

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
            a number between ``0`` and ``1`` is chosen. If ``None``,
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

        Raises
        ------
        :class:`ValueError`
            If the number of `orientations` provided is incorrect.

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
        vertices = []
        edges = []
        choices = [True, False]
        for vertex_id, flip_chance in enumerate(orientations):
            theta = vertex_id*angle_diff
            vertices.append(
                _CycleVertex(
                    id=vertex_id,
                    position=[np.cos(theta), np.sin(theta), 0],
                    flip=generator.choice(
                        choices,
                        p=[flip_chance, 1-flip_chance],
                    ),
                    angle=theta,
                )
            )

            if vertex_id > 0:
                edges.append(
                    Edge(
                        id=len(edges),
                        vertex1=vertices[vertex_id-1],
                        vertex2=vertices[vertex_id],
                    )
                )

        # Save the chosen orientations for __repr__.
        self._orientations = tuple(
            int(vertex.get_flip()) for vertex in vertices
        )

        edges.append(Edge(len(edges), vertices[0], vertices[-1]))
        super().__init__(
            building_block_vertices=self._get_building_block_vertices(
                building_blocks=building_blocks,
                vertices=vertices,
            ),
            edges=tuple(edges),
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=None,
        )

    @staticmethod
    def _normalize_repeating_unit(repeating_unit):
        if isinstance(repeating_unit, tuple):
            return repeating_unit
        base = ord('A')
        return tuple(ord(letter)-base for letter in repeating_unit)

    def _get_building_block_vertices(self, building_blocks, vertices):
        polymer = self._repeating_unit*self._num_repeating_units
        building_block_vertices = {}
        for bb_index, vertex in zip(polymer, vertices):
            bb = building_blocks[bb_index]
            building_block_vertices[bb] = (
                building_block_vertices.get(bb, [])
            )
            building_block_vertices[bb].append(vertex)
        return building_block_vertices

    def _get_scale(self, building_block_vertices):
        length = len(self._repeating_unit)*self._num_repeating_units
        return length*0.25*max(
            bb.get_maximum_diameter()
            for bb in building_block_vertices
        )

    def __repr__(self):
        return (
            f'macrocycle.Macrocycle({self._repeating_unit!r}, '
            f'{self._num_repeating_units!r}, '
            f'{self._orientations!r})'
        )
