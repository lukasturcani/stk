"""
Host Guest Complex
==================

#. :class:`.Complex`

"""

from .vertices import _HostVertex, _GuestVertex
from ..topology_graph import TopologyGraph


class Complex(TopologyGraph):
    """
    Represents a host-guest complex topology graph.

    When using this topology graph, the host must be first in the
    `building_blocks` of the :class:`.ConstructedMolecule`
    and the guest must be second.

    Examples
    --------
    .. code-block:: python

        import stk

        host = stk.ConstructedMolecule(
            building_blocks=(
                stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()]),
                stk.BuildingBlock(
                    smiles='O=CC(C=O)C=O',
                    functional_groups=[stk.AldehydeFactory()],
                ),
            ),
            topology_graph=stk.cage.FourPlusSix()
        )
        host = stk.BuildingBlock.init_from_molecule(host),
        guest = stk.BuildingBlock('[Br][Br]')
        complex1 = stk.ConstructedMolecule(
            building_blocks=(host, guest),
            topology_graph=stk.host_guest.Complex(),
        )

    Change the position and orientation of the guest

    .. code-block:: python

        complex2 = stk.ConstructedMolecule(
            building_blocks=(host, guest),
            topology_graph=stk.host_guest.Complex(
                # Apply a rotation onto the guest molecule such that
                # the vector returned by get_direction() has the same
                # direction as [1, 1, 1].
                guest_start=guest.get_direction(),
                guest_target=[1, 1, 1],
                displacement=[5.3, 2.1, 7.1],
            ),
        )

    """

    def __init__(
        self,
        guest_start=None,
        guest_target=None,
        displacement=(0, 0, 0),
        num_processes=1,
    ):
        """
        Initialize an instance of :class:`.Complex`.

        Parameters
        ----------
        guest_start : :class:`tuple` of :class:`float`, optional
            A direction vector which gets aligned with `guest_target`.

        guest_target : :class:`tuple` of :class:`float`, optional
            A direction vector which determines the rotation applied to
            the guest building block. A rotation such that
            `guest_start` is transformed into `guest_target` is applied
            to the guest building block.

        displacement : :class:`tuple` of :class:`float`, optional
            The translational offset of the guest from the center of
            the host cavity.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        Raises
        ------
        :class:`TypeError`
            If `guest_start` or `guest_target` is defined but the other
            is not.

        """

        num_nones = sum(
            1 for vector in (guest_start, guest_target)
            if vector is None
        )
        if num_nones == 1:
            raise TypeError(
                'If guest_start or guest_target is defined, '
                'the other must be too.'
            )

        if guest_start is None:
            start = target = (1., 0., 0.)
        else:
            start = guest_start = tuple(guest_start)
            target = guest_target = tuple(guest_target)

        # Save the values as None, for __repr__.
        self._guest_start = guest_start
        self._guest_target = guest_target
        self._displacement = displacement

        vertices = (
            _HostVertex(0, [0, 0, 0]),
            _GuestVertex(1, displacement, start, target)
        )
        super().__init__(
            vertices=vertices,
            edges=(),
            reaction_factory=None,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=(),
        )

    def get_building_block_vertices(self, building_blocks):
        return {
            bb: [vertex]
            for bb, vertex in zip(building_blocks, self._vertices)
        }

    def _run_reactions(self, state):
        return state

    def _get_scale(self, building_block_vertices):
        return 1

    def __repr__(self):
        return (
            f'host_guest.Complex('
            f'guest_start={self._guest_start!r}, '
            f'guest_target={self._guest_target!r}, '
            f'displacement={self._displacement!r})'
        )
