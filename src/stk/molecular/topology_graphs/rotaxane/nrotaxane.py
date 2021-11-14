"""
[n]Rotaxane
===========

"""

import numpy as np

from ..topology_graph import NullOptimizer, TopologyGraph
from .vertices import AxleVertex, CycleVertex


class NRotaxane(TopologyGraph):
    """
    Represents [n]rotaxane topology graphs.

    Axle and cycle building blocks do not require functional groups.

    This class assumes one axle with (n-1) macrocycles threaded on it.
    The macrocycles are spaced evenly along the axle in repeating
    patterns. The threaded macrocycles can be described analogously
    to monomers in linear polymers, in terms of a repeating unit,
    except that no bonds are formed between them.

    Examples
    --------
    *Construction*

    .. testcode:: construction

        import stk

        cycle = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='[Br]CC[Br]',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=5,
            ),
        )
        axle = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='AB',
                num_repeating_units=7,
            )
        )
        rotaxane = stk.ConstructedMolecule(
            topology_graph=stk.rotaxane.NRotaxane(
                axle=stk.BuildingBlock.init_from_molecule(axle),
                cycles=(
                    stk.BuildingBlock.init_from_molecule(cycle),
                ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        cycle = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='[Br]CC[Br]',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='A',
                num_repeating_units=5,
            ),
        )
        axle = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='AB',
                num_repeating_units=7,
            )
        )
        rotaxane = stk.ConstructedMolecule(
            topology_graph=stk.rotaxane.NRotaxane(
                axle=stk.BuildingBlock.init_from_molecule(axle),
                cycles=(
                    stk.BuildingBlock.init_from_molecule(cycle),
                ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    rotaxane.get_atoms(),
                    rotaxane.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in rotaxane.get_bonds()
            ),
        )

    *Suggested Optimization*

    For :class:`.NRotaxane` topologies, there is no need to use an
    optimizer, so stick with :class:`.NullOptimizer`. However, it is
    recommended that all building blocks be optimized prior to
    construction.

    *Defining the Orientation of Each Building Block*

    The `orientations` parameter allows the direction of each cycle
    along the axle to be flipped

    .. testcode:: defining-the-orientation-of-each-building-block

        import stk

        cycle1 = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrCNCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=5,
                optimizer=stk.MCHammer(),
            ),
        )
        cycle2 = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrCCNCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=5,
                optimizer=stk.MCHammer(),
            ),
        )
        axle = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        rotaxane1 = stk.ConstructedMolecule(
            topology_graph=stk.rotaxane.NRotaxane(
                axle=stk.BuildingBlock.init_from_molecule(axle),
                cycles=(
                    stk.BuildingBlock.init_from_molecule(cycle1),
                    stk.BuildingBlock.init_from_molecule(cycle2),
                ),
                repeating_unit='AB',
                num_repeating_units=3,
                orientations=(1., 0.5),
            ),
        )

    In the above example, ``cycle1`` is guaranteed to be flipped,
    ``cycle2`` has a 50% chance of being flipped, each time it is
    placed on a node.

    Note that whether a building block will be flipped or not
    is decided during the initialization of :class:`.NRotaxane`

    .. testcode:: defining-the-orientation-of-each-building-block

        # graph1 will always construct the same [n]rotaxane.
        graph1 = stk.rotaxane.NRotaxane(
            axle=stk.BuildingBlock.init_from_molecule(axle),
            cycles=(
                stk.BuildingBlock.init_from_molecule(cycle1),
                stk.BuildingBlock.init_from_molecule(cycle2),
            ),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
        )
        # rotaxane2 and rotaxane3 are guaranteed to be the same as they
        # used the same topology graph.
        rotaxane2 = stk.ConstructedMolecule(graph1)
        rotaxane3 = stk.ConstructedMolecule(graph1)

        # graph2 may lead to a different [n]rotaxane, despite
        # being initialized with the same parameters.
        graph2 = stk.rotaxane.NRotaxane(
            axle=stk.BuildingBlock.init_from_molecule(axle),
            cycles=(
                stk.BuildingBlock.init_from_molecule(cycle1),
                stk.BuildingBlock.init_from_molecule(cycle2),
            ),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45)
        )

        # rotaxane4 and rotaxane5 are guaranteed to be the same because
        # they used the same topology graph. However, they may be
        # different to rotaxane2 and rotaxane3.
        rotaxane4 = stk.ConstructedMolecule(graph2)
        rotaxane5 = stk.ConstructedMolecule(graph2)

    The `random_seed` parameter can be used to get reproducible results

    .. testcode:: defining-the-orientation-of-each-building-block

        # rotaxane6 and rotaxane7 are guaranteed to be the same,
        # because graph3 and graph4 used the same random seed.

        graph3 = stk.rotaxane.NRotaxane(
            axle=stk.BuildingBlock.init_from_molecule(axle),
            cycles=(
                stk.BuildingBlock.init_from_molecule(cycle1),
                stk.BuildingBlock.init_from_molecule(cycle2),
            ),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
            random_seed=4,
        )
        rotaxane6 = stk.ConstructedMolecule(graph3)

        graph4 = stk.rotaxane.NRotaxane(
            axle=stk.BuildingBlock.init_from_molecule(axle),
            cycles=(
                stk.BuildingBlock.init_from_molecule(cycle1),
                stk.BuildingBlock.init_from_molecule(cycle2),
            ),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
            random_seed=4,
        )
        rotaxane7 = stk.ConstructedMolecule(graph4)

    *Using Numbers to Define the Repeating Unit*

    The repeating unit can also be specified through the indices of
    the building blocks

    .. testcode:: using-numbers-to-define-the-repeating-unit

        import stk

        cycle1 = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrCNCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=5,
                optimizer=stk.MCHammer(),
            ),
        )
        cycle2 = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrCCNCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=5,
                optimizer=stk.MCHammer(),
            ),
        )
        cycle3 = stk.ConstructedMolecule(
            topology_graph=stk.macrocycle.Macrocycle(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCNNCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='BrCCNNNCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=5,
                optimizer=stk.MCHammer(),
            ),
        )
        axle = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )

        # r1 and r2 are different ways to write the same thing.
        r1 = stk.ConstructedMolecule(
            topology_graph=stk.rotaxane.NRotaxane(
                axle=stk.BuildingBlock.init_from_molecule(axle),
                cycles=(
                    stk.BuildingBlock.init_from_molecule(cycle1),
                    stk.BuildingBlock.init_from_molecule(cycle2),
                    stk.BuildingBlock.init_from_molecule(cycle3),
                ),
                repeating_unit='ACB',
                num_repeating_units=3,
            ),
        )
        r2 = stk.ConstructedMolecule(
            topology_graph=stk.rotaxane.NRotaxane(
                axle=stk.BuildingBlock.init_from_molecule(axle),
                cycles=(
                    stk.BuildingBlock.init_from_molecule(cycle1),
                    stk.BuildingBlock.init_from_molecule(cycle2),
                    stk.BuildingBlock.init_from_molecule(cycle3),
                ),
                repeating_unit=(0, 2, 1),
                num_repeating_units=3,
            ),
        )

    """

    def __init__(
        self,
        axle,
        cycles,
        repeating_unit,
        num_repeating_units,
        orientations=None,
        random_seed=None,
        num_processes=1,
        optimizer=NullOptimizer(),
    ):
        """
        Initialize a :class:`NRotaxane` instance.

        Parameters
        ----------
        axle : :class:`.BuildingBlock`
            The axle of the rotaxane.

        cycles : :class:`tuple` of :class:`.BuildingBlock`
            The cycles threaded onto the `axle`.

        repeating_unit : :class:`str` or :class:`tuple` of :class:`int`
            A string specifying the repeating unit of the `cycles`.
            For example, ``'AB'`` or ``'ABB'``. The first cycle in
            `cycles` is ``'A'`` and so on.

            The repeating unit can also be specified by the indices of
            `cycles`, for example ``'ABB'`` can be
            written as ``(0, 1, 1)``.

        num_repeating_units : :class:`int`
            The number of repeating units threaded along the axle.

        orientations : :class:`tuple` of :class:`float`, optional
            For each character in the repeating unit, a value
            between ``0`` and ``1`` (both inclusive) must be given in
            a :class:`tuple`. It indicates the probability that each
            cycle will have its orientation along the axle
            flipped. If ``0`` then the cycle is guaranteed not to
            flip. If ``1`` it is guaranteed to flip. This allows the
            user to create head-to-head or head-to-tail chains, as well
            as chain with a preference for head-to-head or head-to-tail
            if a number between ``0`` and ``1`` is chosen. If
            ``None``, then defaults to ``0`` in every case.

            It is also possible to supply an orientation for every
            cycle vertex in the final topology graph. In this case, the
            length of `orientations` must be equal to
            ``len(repeating_unit)*num_repeating_units``.

        random_seed : :class:`int`, optional
            The random seed to use when choosing random orientations.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        optimizer : :class:`.Optimizer`, optional
            Used to optimize the structure of the constructed
            molecule.

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

        vertices = [AxleVertex(0, [0, 0, 0])]
        distance = 1 / (chain_length+1)
        choices = [True, False]
        for vertex_id, flip_chance in enumerate(orientations, 1):
            vertices.append(
                CycleVertex(
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
            optimizer=optimizer,
            edge_groups=None,
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
