"""
Linear
======

"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ....reactions import GenericReactionFactory
from ...topology_graph import (
    Edge,
    NullOptimizer,
    TopologyGraph,
    Vertex,
)
from .vertices import (
    HeadVertex,
    LinearVertex,
    TailVertex,
    UnaligningVertex,
)


class Linear(TopologyGraph):
    """
    Represents a linear polymer topology graph.

    Building blocks with two functional groups are required, unless the
    building block's position is specified to only be at the capping
    positions.

    Examples
    --------
    *Construction*

    Linear polymers require building blocks with two functional groups

    .. testcode:: construction

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=4,
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=4,
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    polymer.get_atoms(),
                    polymer.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in polymer.get_bonds()
            ),
        )


    *Suggested Optimization*

    For :class:`.Linear` topologies, it is recommend to use the
    :class:`.Collapser` optimizer.

    .. testcode:: suggested-optimization

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])

        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=4,
                # Setting scale_steps to False tends to lead to a
                # better structure.
                optimizer=stk.Collapser(scale_steps=False),
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])

        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=4,
                # Setting scale_steps to False tends to lead to a
                # better structure.
                optimizer=stk.Collapser(scale_steps=False),
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    polymer.get_atoms(),
                    polymer.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in polymer.get_bonds()
            ),
        )

    *Construction with Capping Units*

    Building blocks with a single functional group can
    also be provided as capping units

    .. testcode:: construction-with-capping-units

        import stk

        bb1 = stk.BuildingBlock('NCC(F)N', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        bb3 = stk.BuildingBlock('BrCCN', [stk.PrimaryAminoFactory()])

        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2, bb3),
                repeating_unit='ABABC',
                num_repeating_units=1,
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock('NCC(F)N', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        bb3 = stk.BuildingBlock('BrCCN', [stk.PrimaryAminoFactory()])

        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2, bb3),
                repeating_unit='ABABC',
                num_repeating_units=1,
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    polymer.get_atoms(),
                    polymer.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in polymer.get_bonds()
            ),
        )

    *Defining the Orientation of Each Building Block*

    The `orientations` parameter allows the direction of each building
    block along to the chain to be flipped

    .. testcode:: defining-the-orientation-of-each-building-block

        import stk

        bb1 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        bb2 = stk.BuildingBlock(
            smiles='NC(Br)CN',
            functional_groups=[stk.PrimaryAminoFactory()],
        )

        p1 = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=5,
                orientations=(1, 0.5),
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock('O=CC(F)CC=O', [stk.AldehydeFactory()])
        bb2 = stk.BuildingBlock(
            smiles='NC(Br)CN',
            functional_groups=[stk.PrimaryAminoFactory()],
        )

        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=5,
                orientations=(1, 0.5),
                random_seed=1,
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    polymer.get_atoms(),
                    polymer.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in polymer.get_bonds()
            ),
        )

    In the above example, ``bb1`` is guaranteed to be flipped,
    ``bb2`` has a 50% chance of being flipped, each time it is placed
    on a node.

    Note that whether a building block will be flipped or not
    is decided during the initialization of :class:`.Linear`

    .. testcode:: defining-the-orientation-of-each-building-block

        # chain will always construct the same polymer.
        chain = stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
        )
        # p2 and p3 are guaranteed to be the same as they used the same
        # topology graph.
        p2 = stk.ConstructedMolecule(chain)
        p3 = stk.ConstructedMolecule(chain)

        # chain2 may lead to a different polymer than chain, despite
        # being initialized with the same parameters.
        chain2 = stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45)
        )

        # p4 and p5 are guaranteed to be the same because they used
        # the same topology graph. However, they may be different to
        # p2 and p3.
        p4 = stk.ConstructedMolecule(chain2)
        p5 = stk.ConstructedMolecule(chain2)

    The `random_seed` parameter can be used to get reproducible results

    .. testcode:: defining-the-orientation-of-each-building-block

        # p6 and p7 are guaranteed to be the same, because chain3 and
        # chain4 used the same random seed.

        chain3 = stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
            random_seed=4,
        )
        p6 = stk.ConstructedMolecule(chain3)

        chain4 = stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=5,
            orientations=(0.65, 0.45),
            random_seed=4,
        )
        p7 = stk.ConstructedMolecule(chain4)

    *Using Numbers to Define the Repeating Unit*

    The repeating unit can also be specified through the indices of
    the building blocks

    .. testcode:: using-numbers-to-define-the-repeating-unit

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        bb3 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])

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
        optimizer=NullOptimizer(),
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

            If there is only one building block in the constructed
            polymer i.e. the `repeating_unit` has a length of 1 and
            `num_repeating_units` is 1, the building block will not
            be re-oriented, even if you provide a value to
            `orientations`.

        random_seed : :class:`int`, optional
            The random seed to use when choosing random orientations.

        reaction_factory : :class:`.ReactionFactory`, optional
            The factory to use for creating reactions between
            functional groups of building blocks.

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

        polymer_length = len(repeating_unit)*num_repeating_units
        if len(orientations) != polymer_length:
            raise ValueError(
                'The length of orientations must match either '
                'the length of repeating_unit or the '
                'total number of vertices.'
            )

        # Keep these for __repr__.
        self._repeating_unit = self._normalize_repeating_unit(
            repeating_unit=repeating_unit
        )
        self._num_repeating_units = num_repeating_units

        try:
            head, *body, tail = orientations
            vertices_and_edges = self._get_vertices_and_edges(
                head_orientation=head,
                body_orientations=body,
                tail_orientation=tail,
                random_seed=random_seed,
            )
            vertices = vertices_and_edges.vertices
            edges = vertices_and_edges.edges

        except ValueError:
            vertices = (
                UnaligningVertex(0, (0., 0., 0.), False),
            )
            edges = ()

        # Save the chosen orientations for __repr__.
        self._orientations = tuple(int(v.get_flip()) for v in vertices)

        super().__init__(
            building_block_vertices=self._get_building_block_vertices(
                building_blocks=building_blocks,
                vertices=vertices,
            ),
            edges=edges,
            reaction_factory=reaction_factory,
            construction_stages=(),
            optimizer=optimizer,
            num_processes=num_processes,
        )

    @staticmethod
    def _get_vertices_and_edges(
        head_orientation,
        body_orientations,
        tail_orientation,
        random_seed,
    ):
        """
        Get the vertices and edges of the topology graph.

        Parameters
        ----------
        head_orientation : :class:`float`
            The probability that the head vertex will do flipping.

        body_orientations : :class:`iterable` of :class:`float`
            For each body vertex, the probability that it will do
            flipping.

        tail_orientation : :class:`float`
            The probability that the tail vertex will do flipping.

        random_seed : :class:`int`
            The random seed to use.

        Returns
        -------
        :class:`.VerticesAndEdges`
            The vertices and edges of the topology graph.

        """

        generator = np.random.RandomState(random_seed)

        choices = [True, False]
        vertices = [
            HeadVertex(
                id=0,
                position=np.array([0, 0, 0]),
                flip=generator.choice(
                    a=choices,
                    p=[head_orientation, 1-head_orientation],
                ),
            ),
        ]
        edges = []
        for i, p in enumerate(body_orientations, 1):
            flip = generator.choice(choices, p=[p, 1-p])
            vertices.append(
                LinearVertex(i, np.array([i, 0, 0]), flip),
            )
            edges.append(Edge(len(edges), vertices[i-1], vertices[i]))

        vertices.append(
            TailVertex(
                id=len(vertices),
                position=np.array([len(vertices), 0, 0]),
                flip=generator.choice(
                    a=choices,
                    p=[tail_orientation, 1-tail_orientation],
                ),
            ),
        )

        edges.append(Edge(len(edges), vertices[-2], vertices[-1]))

        return _VerticesAndEdges(
            vertices=tuple(vertices),
            edges=tuple(edges),
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

    def _get_building_block_vertices(
        self,
        building_blocks,
        vertices,
    ):
        polymer = self._repeating_unit*self._num_repeating_units
        building_block_vertices = {}
        for bb_index, vertex in zip(polymer, vertices):
            bb = building_blocks[bb_index]
            building_block_vertices[bb] = (
                building_block_vertices.get(bb, [])
            )
            building_block_vertices[bb].append(vertex)

        building_block_vertices = self._with_unaligning_vertices(
            building_block_vertices=building_block_vertices,
        )
        return building_block_vertices

    @staticmethod
    def _with_unaligning_vertices(
        building_block_vertices,
    ):
        clone = {}
        for building_block, vertices in (
            building_block_vertices.items()
        ):
            # Building blocks with 1 placer, cannot be aligned and
            # must therefore use an UnaligningVertex.
            if building_block.get_num_placers() == 1:
                clone[building_block] = tuple(
                    UnaligningVertex(
                        id=vertex.get_id(),
                        position=vertex.get_position(),
                        flip=vertex.get_flip(),
                    ) for vertex in vertices
                )
            else:
                clone[building_block] = vertices

        return clone

    def _get_scale(
        self,
        building_block_vertices,
    ) -> float:
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


@dataclass(frozen=True)
class _VerticesAndEdges:
    vertices: tuple[Vertex]
    edges: tuple[Edge]
