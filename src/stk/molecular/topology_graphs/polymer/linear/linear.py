"""
Linear
======

"""

import numpy as np

from .vertices import _HeadVertex, _TailVertex, _LinearVertex
from ...topology_graph import TopologyGraph, Edge
from ....reactions import GenericReactionFactory


class Linear(TopologyGraph):
    """
    Represents a linear polymer topology graph.

    Examples
    --------
    *Construction*

    Linear polymers require building blocks with two functional groups

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=12,
            ),
        )

    *Construction with Capping Units*

    Building blocks with a single functional group can
    also be provided as capping units

    .. code-block:: python

        bb3 = stk.BuildingBlock('CCN', [stk.PrimaryAminoFactory()])
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2, bb3),
                repeating_unit='ABABC',
                num_repeating_units=1,
            ),
        )

    *Defining the Orientation of Each Building Block*

    The `orientations` parameter allows the direction of each building
    block along to the chain to be flipped

    .. code-block:: python

        bb4 = stk.BuildingBlock('NCNCCN', [stk.PrimaryAminoFactory()])

        p3 = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb2, bb4),
                repeating_unit='AB',
                num_repeating_units=5,
                orientations=(1, 0.5),
            ),
        )

    In the above example, ``bb2`` is guaranteed to be flipped,
    ``bb4`` has a 50% chance of being flipped, each time it is placed
    on a node.

    Note that whether a building block will be flipped or not
    is decided during the initialization of :class:`.Linear`

    .. code-block:: python

        # chain will always construct the same polymer.
        chain = stk.polymer.Linear(
            building_blocks=(bb2, bb4),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
        )
        # p4 and p5 are guaranteed to be the same as they used the same
        # topology graph.
        p4 = stk.ConstructedMolecule(chain)
        p5 = stk.ConstructedMolecule(chain)

        # chain2 may lead to a different polymer than chain, despite
        # being initialized with the same parameters.
        chain2 = stk.polymer.Linear(
            building_blocks=(bb2, bb4),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45)
        )

        # p6 and p7 are guaranteed to be the same because they used
        # the same topology graph. However, they may be different to
        # p4 and p5.
        p6 = stk.ConstructedMolecule(chain2)
        p7 = stk.ConstructedMolecule(chain2)

    The `random_seed` parameter can be used to get reproducible results

    .. code-block:: python

        # p8 and p9 are guaranteed to be the same, because chain3 and
        # chain4 used the same random seed.

        chain3 = stk.polymer.Linear(
            building_blocks=(bb2, bb4),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
            random_seed=4,
        )
        p8 = stk.ConstructedMolecule(chain3)

        chain4 = stk.polymer.Linear(
            building_blocks=(bb2, bb4),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
            random_seed=4,
        )
        p9 = stk.ConstructedMolecule(chain4)

    *Using Numbers to Define the Repeating Unit*

    The repeating unit can also be specified through the indices of
    the building blocks

    .. code-block:: python

        # p1 and p2 are different ways to write the same thing.
        p1 = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2, bb3),
                repeating_unit='ACB',
                num_repeating_units=1,
            ),
        )
        p2 = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2, bb3),
                repeating_unit=(0, 2, 1),
                num_repeating_units=1,
            ),
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
        num_processes=1,
    ):
        """
        Initialize a :class:`Linear` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks of the polymer.

        repeating_unit : :class:`str` or :class:`tuple` of :class:`int`
            A string specifying the repeating unit of the polymer.
            For example, ``'AB'`` or ``'ABB'``. The first building
            block passed to `building_blocks` is ``'A'`` and so on.

            The repeating unit can also be specified by the indices of
            `building_blocks`, for example ``'ABB'`` can be
            written as ``(0, 1, 1)``.

        num_repeating_units : :class:`int`
            The number of repeating units which are used to make the
            polymer.

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
            then ``0`` is picked in all cases.

            It is also possible to supply an orientation for every
            vertex in the final topology graph. In this case, the
            length of `orientations` must be equal to
            ``len(repeating_unit)*num_repeating_units``.

        random_seed : :class:`int`, optional
            The random seed to use when choosing random orientations.

        reaction_factory : :class:`.ReactionFactory`, optional
            The factory to use for creating reactions between
            functional groups of building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        Raises
        ------
        :class:`ValueError`
            If the length of `orientations` is not equal in length to
            `repeating_unit` or to the total number of vertices.

        """

        if orientations is None:
            orientations = tuple(
                0. for i in range(len(repeating_unit))
            )

        if len(orientations) == len(repeating_unit):
            orientations = orientations*num_repeating_units

        polymer_length = len(repeating_unit)*num_repeating_units
        if len(orientations) != polymer_length:
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

        head, *body, tail = orientations
        choices = [True, False]
        vertices = [
            _HeadVertex(
                id=0,
                position=np.array([0, 0, 0]),
                flip=generator.choice(choices, p=[head, 1-head])
            ),
        ]
        edges = []
        for i, p in enumerate(body, 1):
            flip = generator.choice(choices, p=[p, 1-p])
            vertices.append(
                _LinearVertex(i, np.array([i, 0, 0]), flip)
            )
            edges.append(Edge(len(edges), vertices[i-1], vertices[i]))

        vertices.append(
            _TailVertex(
                id=len(vertices),
                position=np.array([len(vertices), 0, 0]),
                flip=generator.choice(choices, p=[tail, 1-tail])),
        )

        # Save the chosen orientations for __repr__.
        self._orientations = tuple(int(v.get_flip()) for v in vertices)

        edges.append(Edge(len(edges), vertices[-2], vertices[-1]))

        super().__init__(
            building_block_vertices=self._get_building_block_vertices(
                building_blocks=building_blocks,
                vertices=vertices,
            ),
            edges=tuple(edges),
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
        )

    def clone(self):
        clone = super().clone()
        clone._repeating_unit = self._repeating_unit
        clone._num_repeating_units = self._num_repeating_units
        clone._orientations = self._orientations
        return clone

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
        return max(
            bb.get_maximum_diameter()
            for bb in building_block_vertices
        )

    def __repr__(self):
        return (
            f'polymer.Linear({self._repeating_unit!r}, '
            f'{self._num_repeating_units!r}, '
            f'{self._orientations!r})'
        )
