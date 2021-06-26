"""
Host Guest Complex
==================

"""

from ...molecules import BuildingBlock
from .vertices import HostVertex, GuestVertex
from ..topology_graph import TopologyGraph, NullOptimizer
from dataclasses import dataclass


@dataclass(frozen=True)
class Guest:
    building_block: BuildingBlock
    start_vector: tuple = (1., 0., 0.)
    end_vector: tuple = (1., 0., 0.)
    displacement: tuple = (0., 0., 0.)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'guest_start={self.start_vector!r}, '
            f'guest_target={self.end_vector!r}, '
            f'displacement={self.displacement!r})'
        )


class Complex(TopologyGraph):
    """
    Represents a host-guest complex topology graph.

    Host and guest building blocks do not require functional groups.

    Examples
    --------
    *Construction*

    You can use :class:`.ConstructedMolecule` instances as the host,
    but you should turn them into a :class:`.BuildingBlock` first

    .. testcode:: construction

        import stk

        host = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='BrCC(Br)CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
            ),
        )
        complex = stk.ConstructedMolecule(
            topology_graph=stk.host_guest.Complex(
                host=stk.BuildingBlock.init_from_molecule(host),
                guest=(stk.host_guest.Guest(
                    building_block=stk.BuildingBlock('[Br][Br]'),
                ), ),
            ),
        )

    *Suggested Optimization*

    For :class:`.Complex` topologies, there is no need to use an
    optimizer, so stick with :class:`.NullOptimizer`. However, it is
    recommended that all building blocks be optimized prior to
    construction.

    *Changing the Position of the Guest*

    You can change the position and orientation of the guest, as well
    as its displacement

    .. testcode:: changing-the-position-of-the-guest

        import stk

        host = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='BrCC(Br)CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
            ),
        )
        guest = stk.host_guest.Guest(stk.BuildingBlock('[Br][Br]'))
        complex = stk.ConstructedMolecule(
            topology_graph=stk.host_guest.Complex(
                host=stk.BuildingBlock.init_from_molecule(host),
                guests=(guest, ),
                # Apply a rotation onto the guest molecule such that
                # the vector returned by get_direction() has the same
                # direction as [1, 1, 1].
                guest_start=guest.get_direction(),
                guest_target=[1, 1, 1],
                # Change the displacement of the guest.
                displacement=[5.3, 2.1, 7.1],
            ),
        )

    """

    def __init__(
        self,
        host,
        guests,
        num_processes=1,
        optimizer=NullOptimizer(),
    ):
        """
        Initialize an instance of :class:`.Complex`.

        Parameters
        ----------
        host : :class:`.BuildingBlock`
            The host molecule.

        guest : :class:`.BuildingBlock`
            The guest molecule.

        guest_start : :class:`tuple` of :class:`float`, optional
            A direction vector which gets aligned with `guest_target`.

        guest_target : :class:`tuple` of :class:`float`, optional
            A direction vector which determines the rotation applied to
            the `guest`. A rotation such that `guest_start` is
            transformed into `guest_target` is applied.

        displacement : :class:`tuple` of :class:`float`, optional
            The translational offset of the guest.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        optimizer : :class:`.Optimizer`, optional
            Used to optimize the structure of the constructed
            molecule.

        Raises
        ------
        :class:`TypeError`
            If `guest_start` or `guest_target` is defined but the other
            is not.

        """

        building_block_vertices = {host: (HostVertex(0, [0, 0, 0]), )}
        guest_vertices = {
            guest.building_block: (GuestVertex(
                id=i+1,
                position=guest.displacement,
                start=guest.start_vector,
                target=guest.end_vector,
            ), )
            for i, guest in enumerate(guests)
        }
        building_block_vertices.update(guest_vertices)

        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=(),
            reaction_factory=None,
            construction_stages=(),
            num_processes=num_processes,
            optimizer=optimizer,
            edge_groups=(),
        )

    def clone(self):
        clone = super().clone()
        return clone

    def _run_reactions(self, state):
        return state

    def _get_scale(self, building_block_vertices):
        return 1

    def __repr__(self):
        return (
            'host_guest.Complex('
            # f'guest_start={self._guest_start!r}, '
            # f'guest_target={self._guest_target!r}, '
            # f'displacement={self._displacement!r})'
        )
