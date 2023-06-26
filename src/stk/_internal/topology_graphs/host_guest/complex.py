from __future__ import annotations

import typing
from collections import abc

from stk._internal.building_block import BuildingBlock
from stk._internal.construction_state import (
    ConstructionState,
)
from stk._internal.optimizers.null import (
    NullOptimizer,
)
from stk._internal.optimizers.optimizer import (
    Optimizer,
)
from stk._internal.reaction_factories.generic_reaction_factory import (
    GenericReactionFactory,
)
from stk._internal.topology_graphs.topology_graph.topology_graph import (
    TopologyGraph,
)
from stk._internal.topology_graphs.vertex import (
    Vertex,
)

from .vertices import GuestVertex, HostVertex


class Guest:
    """
    Holds the data defining the placement of a guest molecule.

    """

    def __init__(
        self,
        building_block: BuildingBlock,
        start_vector: tuple[float, float, float] = (1.0, 0.0, 0.0),
        end_vector: tuple[float, float, float] = (1.0, 0.0, 0.0),
        displacement: tuple[float, float, float] = (1.0, 0.0, 0.0),
    ) -> None:
        """
        Initialize a :class:`.Guest` instance.

        Parameters:
            building_block: The guest molecule.

            start_vector: A direction vector which gets aligned with
                `end_vector`.

            end_vector: A direction vector which determines the
                rotation applied to the `building_block`. A rotation
                such that `start_vector` is transformed into
                `end_vector` is applied.

            displacement: The translational offset of the guest.

        """

        self._building_block = building_block
        self._start_vector = start_vector
        self._end_vector = end_vector
        self._displacement = displacement

    def get_building_block(self) -> BuildingBlock:
        """
        Return the building block.

        Returns:
            The building block.

        """

        return self._building_block

    def get_start_vector(self) -> tuple[float, float, float]:
        """
        Return the start vector.

        Returns:
            The start vector.

        """

        return self._start_vector

    def get_end_vector(self) -> tuple[float, float, float]:
        """
        Return the end vector.

        Returns:
            The end vector.

        """

        return self._end_vector

    def get_displacement(self) -> tuple[float, float, float]:
        """
        Return the displacement.

        Returns:
            The displacement.

        """

        return self._displacement

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"{self._building_block!r}, "
            f"start_vector={self._start_vector!r}, "
            f"end_vector={self._end_vector!r}, "
            f"displacement={self._displacement!r})"
        )


class Complex(TopologyGraph):
    """
    Represents a host-guest complex topology graph.

    Host and guest building blocks do not require functional groups.

    Examples:

        *Construction*

        You can use :class:`.ConstructedMolecule` instances as the
        host, but you should turn them into a :class:`.BuildingBlock`
        first

        .. testcode:: construction

            import stk

            host = stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='NC1CCCCC1N',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                        stk.BuildingBlock(
                            smiles='O=Cc1cc(C=O)cc(C=O)c1',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    ),
                    optimizer=stk.MCHammer(),
                ),
            )
            complex = stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock.init_from_molecule(host),
                    guests=stk.host_guest.Guest(
                        building_block=stk.BuildingBlock('[Br][Br]'),
                    ),
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            host = stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='NC1CCCCC1N',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                        stk.BuildingBlock(
                            smiles='O=Cc1cc(C=O)cc(C=O)c1',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    ),
                    optimizer=stk.MCHammer(),
                ),
            )
            complex = stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock.init_from_molecule(host),
                    guests=stk.host_guest.Guest(
                        building_block=stk.BuildingBlock('[Br][Br]'),
                    ),
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        complex.get_atoms(),
                        complex.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in complex.get_bonds()
                ),
            )

        You can also generate complexes with multiple guests.

        .. testcode:: multi-guest-construction

            import stk

            host = stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='NC1CCCCC1N',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                        stk.BuildingBlock(
                            smiles='O=Cc1cc(C=O)cc(C=O)c1',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    ),
                    optimizer=stk.MCHammer(),
                ),
            )
            guest1 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('BrBr'),
                displacement=(0., 3., 0.),
            )
            guest2 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('C1CCCC1'),
            )

            complex = stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock.init_from_molecule(host),
                    guests=(guest1, guest2),
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            host = stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='NC1CCCCC1N',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                        stk.BuildingBlock(
                            smiles='O=Cc1cc(C=O)cc(C=O)c1',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    ),
                    optimizer=stk.MCHammer(),
                ),
            )
            guest1 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('BrBr'),
                displacement=(0., 3., 0.),
            )
            guest2 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('C1CCCC1'),
            )

            complex = stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock.init_from_molecule(host),
                    guests=(guest1, guest2),
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        complex.get_atoms(),
                        complex.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in complex.get_bonds()
                ),
            )

        *Suggested Optimization*

        For :class:`.Complex` topologies, it is recommended to use the
        :class:`.Spinner` optimizer. It is also recommended that the
        building blocks are already optimized prior to construction.
        This optimizer will work on multi-guest systems.

        .. testcode:: suggested-optimization

            import stk

            host = stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='NC1CCCCC1N',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                        stk.BuildingBlock(
                            smiles='O=Cc1cc(C=O)cc(C=O)c1',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    ),
                    optimizer=stk.MCHammer(),
                ),
            )
            guest1 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('BrBr'),
                displacement=(0., 3., 0.),
            )
            guest2 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('C1CCCC1'),
            )

            complex = stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock.init_from_molecule(host),
                    guests=(guest1, guest2),
                    optimizer=stk.Spinner(),
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            host = stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='NC1CCCCC1N',
                            functional_groups=[
                                stk.PrimaryAminoFactory(),
                            ],
                        ),
                        stk.BuildingBlock(
                            smiles='O=Cc1cc(C=O)cc(C=O)c1',
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    ),
                    optimizer=stk.MCHammer(),
                ),
            )
            guest1 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('BrBr'),
                displacement=(0., 3., 0.),
            )
            guest2 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('C1CCCC1'),
            )

            complex = stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock.init_from_molecule(host),
                    guests=(guest1, guest2),
                    optimizer=stk.Spinner(),
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        complex.get_atoms(),
                        complex.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in complex.get_bonds()
                ),
            )

        *Changing the Position of the Guest*

        You can change the position and orientation of the guest, as
        well as its displacement

        .. testcode:: changing-the-position-of-the-guest

            import stk

            host = stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCC(Br)CBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                ),
            )

            guest_building_block = stk.BuildingBlock('[Br][Br]')
            guest = stk.host_guest.Guest(
                building_block=guest_building_block,
                # Apply a rotation onto the guest molecule such that
                # the vector returned by get_direction() has the same
                # direction as [1, 1, 1].
                start_vector=guest_building_block.get_direction(),
                end_vector=[1, 1, 1],
                # Change the displacement of the guest.
                displacement=[5.3, 2.1, 7.1],
            )
            complex = stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock.init_from_molecule(host),
                    guests=guest,
                ),
            )

    """

    def __init__(
        self,
        host: BuildingBlock,
        guests: typing.Union[Guest, typing.Iterable[Guest]],
        num_processes: int = 1,
        optimizer: Optimizer = NullOptimizer(),
    ) -> None:
        """
        Initialize an instance of :class:`.Complex`.

        Parameters:
            host: The host molecule.

            guests: The guest molecules. Can be a single
                :class:`.Guest` instance if only one guest is being
                used.

            num_processes: The number of parallel processes to create
                during :meth:`construct`.

            optimizer: Used to optimize the structure of the
                constructed molecule.

        """

        building_block_vertices = self._get_vertices_from_guests(
            host=host,
            guests=guests,
        )

        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=(),
            reaction_factory=GenericReactionFactory(),
            construction_stages=(),
            num_processes=num_processes,
            optimizer=optimizer,
            edge_groups=(),
        )

    def _get_vertices_from_guests(
        self,
        host: BuildingBlock,
        guests: typing.Union[Guest, typing.Iterable[Guest]],
    ) -> dict[BuildingBlock, abc.Sequence[Vertex]]:
        if isinstance(guests, Guest):
            guests = (guests,)

        building_block_vertices: dict[BuildingBlock, abc.Sequence[Vertex]]
        building_block_vertices = {host: (HostVertex(0, (0.0, 0.0, 0.0)),)}
        guest_vertices = {
            guest.get_building_block(): (
                GuestVertex(
                    id=i + 1,
                    position=guest.get_displacement(),
                    start=guest.get_start_vector(),
                    target=guest.get_end_vector(),
                ),
            )
            for i, guest in enumerate(guests)
        }
        building_block_vertices.update(guest_vertices)

        return building_block_vertices

    def clone(self) -> Complex:
        return self._clone()

    def _run_reactions(
        self,
        state: ConstructionState,
    ) -> ConstructionState:
        return state

    def _get_scale(
        self,
        building_block_vertices: dict[BuildingBlock, abc.Sequence[Vertex]],
    ) -> float:
        return 1.0

    def __repr__(self) -> str:
        return "host_guest.Complex()"
