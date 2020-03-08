import numpy as np

from .vertices import _AxleVertex, _CycleVertex
from ..topology_graph import TopologyGraph


class NRotaxane(TopologyGraph):
    """
    Represents [n]rotaxane topology graphs.

    This class assumes one axle with (n-1) macrocycles threaded on it.
    The macrocycles are spaced evenly along the thread in repeating
    patterns. The threaded macrocycles can be described analogously
    to monomers in linear polymers, in terms of a repeating unit,
    except that no bonds are formed between them.

    The axle must be provided first to the `building_blocks` in
    :class:`.ConstructedMolecule.__init__`.

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

    :class:`.NRotaxane` shares many parameters with :class:`.Linear`,
    and the examples described there are also valid for this class.
    Be sure to read them.

    """

    def __init__(
        self,
        axle,
        cycles,
        repeating_unit,
        num_repeating_units,
        orientations=None,
        random_seed=None,
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

            It is also possible to supply an orientation for every
            cycle vertex in the final topology graph. In this case, the
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

        self._repeating_unit = self._normalize_repeating_unit(
            repeating_unit=repeating_unit,
        )
        self._num_repeating_units = num_repeating_units

        vertices = [_AxleVertex(0, [0, 0, 0])]
        distance = 1 / (chain_length+1)
        choices = [True, False]
        for vertex_id, flip_chance in enumerate(orientations, 1):
            vertices.append(
                _CycleVertex(
                    id=vertex_id,
                    position=[vertex_id*distance-0.5, 0, 0],
                    flip=generator.choice(
                        choices,
                        p=[flip_chance, 1-flip_chance],
                    )
                )
            )

        # Save the chosen orientations for __repr__.
        self._orientations = tuple(
            int(vertex.get_flip()) for vertex in vertices[1:]
        )

        super().__init__(
            building_block_vertices=self._get_building_block_vertices(
                axle=axle,
                cycles=cycles,
                vertices=vertices,
            ),
            edges=(),
            reaction_factory=None,
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

    def _get_building_block_vertices(self, axle, cycles, vertices):
        threads = self._repeating_unit*self._num_repeating_units
        building_block_vertices = {}
        building_block_vertices[axle] = vertices[0:1]
        for cycle_index, vertex in zip(threads, vertices[1:]):
            bb = cycles[cycle_index]
            building_block_vertices[bb] = (
                building_block_vertices.get(bb, [])
            )
            building_block_vertices[bb].append(vertex)
        return building_block_vertices

    def _run_reactions(self, state):
        return state

    def _get_scale(self, building_block_vertices):
        axle = next(iter(building_block_vertices))
        return 0.8*axle.get_maximum_diameter()

    def __repr__(self):
        return (
            f'rotaxane.NRotaxane('
            f'{self._repeating_unit!r}, '
            f'{self._num_repeating_units}, '
            f'{self._orientations!r}'
            f')'
        )
